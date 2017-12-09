# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 19:59:10 2017

@author: Alpha
"""

import numpy as np
import scipy as sc
from sklearn.linear_model import LogisticRegression as LR
from sklearn.model_selection import train_test_split
from sklearn import svm

lin_clf = svm.LinearSVC()
#lin_clf = svm.SVC()
lr = LR()
data = sc.io.loadmat('E:/data/HippoWaveleteAD_NC_811_412.mat')
#data = sc.io.loadmat('E:/data/HippoWaveleteAD_MCI_811_585.mat')

X = data['D']
X = X.T
Y = data['FD']

pres = []
for j in range(100):
    pre = []
    for i in range(7):
        if i == 6 :
            X2 = X
        else:
            X2 = X[:,26*i:(i+1)*26]
        
        
        X_train, X_test, y_train, y_test = train_test_split(X2, Y,
                                                            random_state=43,
                                                            stratify=Y,
                                                            test_size=0.97)
        
        lr.fit(X_train,y_train)
        lin_clf.fit(X_train,y_train)
        
        #print(lr.score(X_test,y_test))
        pre.append(lin_clf.score(X_test,y_test))
    pres.append(pre)
#        print(lin_clf.score(X_test,y_test))
pres = np.array(pres)
100*pres.mean(axis = 0)
#prediction = prediction.astype(np.float) 
#y_test = y_test.astype(np.float)
#y_test = y_test.reshape(400)