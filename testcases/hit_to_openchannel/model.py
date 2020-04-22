#!/usr/bin/env python2

import numpy as np
from sklearn.linear_model import LinearRegression, Lasso, LogisticRegression
from sklearn.cross_validation import train_test_split
import sys

np.set_printoptions(linewidth=np.nan)

# =============================================================================
# Linear regression
# =============================================================================

# def predict(X, y):
#     # X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=0.2)
#     X_train, X_test, y_train, y_test = X, X, y, y
#     for alpha in [1000, 100, 10, 1, 0.1, 0.01, 0.001, 0.0001, 1e-5, 1e-6]:
#         model = Lasso(alpha=alpha, normalize=True).fit(X_train, y_train)
#         score = model.score(X_test, y_test)
#         print 'alpha = %f\tR^2 = %.10f\tcoefs = %s' % (alpha, score, model.coef_)
#     model = LinearRegression(normalize=True).fit(X_train, y_train)
#     score = model.score(X_test, y_test)
#     print 'alpha = %f\tR^2 = %.10f\tcoefs = %s' % (0.0, score, model.coef_)

# Algebraic & DOM separate:
# ALL = np.genfromtxt('all.csv', dtype='S20,<i8,<i8,<i8,<i8,<i8,<i8,<i8,<i8,?,<i8,<i8,<f8,<f8', delimiter=',', names=True)
# A = ALL[ np.logical_not(ALL['use_dom']) ]
# X = np.vstack((A['flow_x'], A['flow_y'], A['flow_z'], A['parcel_size'], A['rk_order'], A['ftts'])).T
# print 'Algebraic rho:'
# predict(X, A['rho'])
# print 'Algebraic cost:'
# predict(X, A['cost'])
# A = ALL[ ALL['use_dom'] ]
# X = np.vstack((A['flow_x'], A['flow_y'], A['flow_z'], A['parcel_size'], A['dom_x'], A['dom_y'], A['dom_z'], A['quads'], A['rk_order'], A['ftts'])).T
# print 'DOM rho:'
# predict(X, A['rho'])
# print 'DOM cost:'
# predict(X, A['cost'])

# Algebraic & DOM together:
# A = np.genfromtxt('all2.csv', dtype='S20,<i8,<i8,<i8,<i8,<i8,<i8,<i8,<i8,<i8,<i8,<f8,<f8', delimiter=',', names=True)
# X = np.vstack((A['flow_x'], A['flow_y'], A['flow_z'], A['parcel_size'], A['dom_x'], A['dom_y'], A['dom_z'], A['quads'], A['rk_order'], A['ftts'])).T
# print 'rho:'
# predict(X, A['rho'])
# print 'cost:'
# predict(X, A['cost'])

# =============================================================================
# Logistic regression
# =============================================================================

class Result(object):
    def __init__(self, X, r, threshold, train_size):
        y = r > threshold
        X_train, X_test, r_train, r_test = train_test_split(X, r, train_size=train_size)
        y_train = r_train > threshold
        y_test = r_test > threshold
        if sum(y_train) == len(y_train):
            y_pred = np.array([True] * len(y_test))
        elif sum(y_train) == 0:
            y_pred = np.array([False] * len(y_test))
        else:
            model = LogisticRegression().fit(X_train, y_train)
            y_pred = model.predict(X_test)
        self.threshold = threshold
        self.train_size = train_size
        self.pop_pos_rate = float(sum(y)) / len(y)
        self.pred_pos_rate = float(sum(y_pred)) / len(y_pred)
        self.work_saved = 1.0 - self.pred_pos_rate*(1.0-train_size) - train_size
        self.misc_pos = sum([real and not pred for (real,pred) in zip(y_test,y_pred)])
        self.real_pos = sum(y_test)
        self.fn_rate = float(self.misc_pos)/self.real_pos if self.real_pos != 0 else 0.0
        self.misc_neg = sum([not real and pred for (real,pred) in zip(y_test,y_pred)])
        self.real_neg = sum(np.logical_not(y_test))
        self.fp_rate = float(self.misc_neg)/self.real_neg if self.real_neg != 0 else 0.0
        self.corr_loss = max(r) - max(np.concatenate((r_train, r_test[y_pred.nonzero()])))

    @staticmethod
    def print_header(self):
        print('Response threshold' + '\t' +
              'Training set size' + '\t' +
              'Population acceptance rate' + '\t' +
              'Classifier acceptance rate' + '\t' +
              'Work saved' + '\t' +
              'FNs' + '\t' +
              'FN rate' + '\t' +
              'FPs' + '\t' +
              'FP rate' + '\t' +
              'Correlation loss')

    def print_self(self):
        print('{:.3f}'.format(self.threshold) + '\t' +
              '{:.2f}%'.format(100.0*self.train_size) + '\t' +
              '{:.2f}%'.format(100.0*self.pop_pos_rate) + '\t' +
              '{:.2f}%'.format(100.0*self.pred_pos_rate) + '\t' +
              '{:.2f}%'.format(100.0*self.work_saved) + '\t' +
              '"{:d}/{:d}"'.format(self.misc_pos, self.real_pos) + '\t' +
              '{:.2f}%'.format(100.0*self.fn_rate) + '\t' +
              '"{:d}/{:d}"'.format(self.misc_neg, self.real_neg) + '\t' +
              '{:.2f}%'.format(100.0*self.fp_rate) + '\t' +
              '{:.3f}'.format(self.corr_loss))

A = np.genfromtxt('all2.csv', dtype='S20,<i8,<i8,<i8,<i8,<i8,<i8,<i8,<i8,<i8,<i8,<f8,<f8', delimiter=',', names=True)
X = np.vstack((A['flow_x'], A['flow_y'], A['flow_z'], A['parcel_size'], A['dom_x'], A['dom_y'], A['dom_z'], A['quads'], A['rk_order'], A['ftts'])).T
r = A['rho']
print 'threshold train_size work_saved corr_loss'
for threshold in [x/200.0 for x in range(150,200)]:
    for train_size in [x/20.0 for x in range(1,20)]:
        results = [Result(X, r, threshold, train_size) for dummy in range(0,100)]
        print threshold, train_size, np.mean([res.work_saved for res in results]), np.mean([res.corr_loss for res in results])
