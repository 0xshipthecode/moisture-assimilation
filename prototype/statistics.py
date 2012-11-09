
import numpy as np


def compute_ols_estimator(X, y):
    """
    Compute the ordinary least squares estimator of y using x.
    X has variables in columns, observations in rows.  Y is 1D array of observations.
    
      beta = compute_ols_estimator(x, y)
      
    """
    Nb = X.shape[1]
    Xm = np.asmatrix(X)
    Ym = np.asmatrix(y)
    beta = np.linalg.inv(Xm.T * Xm) * Xm.T * Ym
    return np.asarray(beta)[:,0]
