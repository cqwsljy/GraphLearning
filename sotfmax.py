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
clf = svm.SVC()
lr = LR()
data = sc.io.loadmat('E:/data/HippoWaveleteAD_NC_811_412.mat')
X = data['D']
X = X.T
Y = data['FD']

X_train, X_test, y_train, y_test = train_test_split(X, Y,
                                                    random_state=42,
                                                    stratify=Y,
                                                    test_size=0.2)

lr.fit(X_train,y_train)
lr.score(X_test,y_test)

clf.fit(X_train,y_train)
clf.score(X_test,y_test)