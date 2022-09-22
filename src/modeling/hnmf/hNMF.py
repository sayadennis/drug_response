import numpy as np
import pandas as pd
import torch
from torch import nn
from collections import defaultdict
from sklearn.metrics import accuracy_score, balanced_accuracy_score
import utils

## assumes X is passed as [Xtrain; Xtest], and only y_train is passed
class hNMF(nn.Module):
    def __init__(self, X1, X2, 
            k = 10, n_iter = 10, eps = 1e-7, whybrid = 1.,
            floss1 = 'l2', floss2 = 'kl', weight_decay = 1e-5, wortho=1, C=1,
            lr = 1e-2, verbose = False, seed=None, fn=None,
            device=torch.device('cpu')):
        super(hNMF, self).__init__()
        if seed is not None:
            torch.manual_seed(seed)
            np.random.seed(seed)
        self.n_iter = n_iter
        self.k = k
        self.floss1 = floss1
        self.floss2 = floss2
        self.weight_decay = weight_decay
        self.lr = lr
        self.verbose = verbose
        self.eps = eps
        self.wortho = wortho
        self.fn = fn
        self.C = C
        self.decomposed = False
        self.whybrid = 1. # tradeoff weight of the hybrid NMF between the two matrices
        self.device = device
        self.report = defaultdict(list)
        self.__initfact__(X1, X2)

    def __initfact__(self, X1, X2):
        # Define data dimensions as class attr
        self.n = torch.tensor(X1.shape, dtype=int)[0].to(self.device)
        self.m1 = torch.tensor(X1.shape, dtype=int)[1].to(self.device)
        self.m2 = torch.tensor(X2.shape, dtype=int)[1].to(self.device)

        # Define input matrices as class attr
        self.X1 = torch.from_numpy(X1).float().to(self.device)
        self.X2 = torch.from_numpy(X2).float().to(self.device)
        
        # Random initialization of the 
        self.scale = torch.tensor(torch.mean(torch.cat((self.X1, self.X2), dim=1)) / self.k).to(self.device)
        W = torch.abs(torch.rand([self.n,self.k]).to(self.device) * self.scale).to(self.device)
        H1 = torch.abs(torch.rand([self.k,self.m1])).to(self.device)
        H2 = torch.abs(torch.rand([self.k,self.m2])).to(self.device)
        # self.scale = torch.tensor(torch.mean(self.X) / self.k).to(self.device)

        # Define matrices as torch.nn.Parameter() 
        self.W = torch.nn.Parameter(W)
        self.H1 = torch.nn.Parameter(H1)
        self.H2 = torch.nn.Parameter(H2)
        self.identity = torch.eye(self.k, device=self.device)
        
        # set factorization loss 1 
        if self.floss1 == 'l2':
            self.loss_fac1 = utils.l2
        elif self.floss1 == 'kl':
            self.loss_fac1 = utils.kl_div
        # set factorization loss 2 
        if self.floss2 == 'l2':
            self.loss_fac2 = utils.l2
        elif self.floss2 == 'kl':
            self.loss_fac2 = utils.kl_div

        # define optimizer 
        self.opt = torch.optim.Adam(self.parameters(), lr=self.lr, weight_decay=self.weight_decay)

    def to(self, device):
        self.device = device
        self.X1 = self.X1.to(device)
        self.X2 = self.X2.to(device)
        return super(hNMF, self).to(device)

    def plus(self, X):
        X[X < 0] = 0 # self.eps
        return X

    def __autograd__(self, epoch):
        """
           autograd update, with gradient projection
        """
        self.opt.zero_grad()
        l = self.loss_fac1(self.W @ self.H1, self.X1) + self.whybrid * self.loss_fac2(self.W @ self.H2, self.X2) # * self.k * self.n
        #### ADD LINE TO COMBINE LOSS 1 AND 2 #### 
        # add l2 regularization for orthogonality
        WtW = torch.mm(torch.t(self.W), self.W)
        if torch.mean(WtW) > self.eps:
            WtW = WtW/torch.mean(WtW)
        ## THE BELOW MIGHT NEED EDITING ## - using L2 loss for this but not sure if this is appropriate 
        l += utils.l2(WtW/self.k, self.identity) * self.wortho * self.k # scale invariant orthogonal

        l.backward()
        self.opt.step()
        ## grad projection
        self.W.data = self.plus(self.W.data)
        self.H1.data = self.plus(self.H1.data)
        self.H2.data = self.plus(self.H2.data)
        return l.item()

    def fit(self):
        for e in range(self.n_iter):
            l = self.__autograd__(e)

            self.report['epoch'].append(e)
            self.report['loss'].append(l)
            if self.verbose and e % 100 == 0:
                print("%d\tloss: %.4f"%(e,l))
        self.decomposed = True
        return self

    def show_report(self):
        return pd.DataFrame(self.report)

    def fit_transform(self):
        if not self.decomposed:
            self.fit()
            # detach all params including the linear fc layer
            for p in self.parameters():
                p.requires_grad = False
        if self.device.type == 'cuda':
            return [self.W.detach().cpu().numpy(), self.H1.detach().cpu().numpy(), self.H2.detach().cpu().numpy()]
        else:
            return [self.W.detach().numpy(), self.H1.detach().numpy(), self.H2.detach().numpy()]
