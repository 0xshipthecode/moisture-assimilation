
import numpy as np
from diagnostics import diagnostics


class MeanFieldModel:
    """
    The model of the mean field for simple kriging.  The idea is to take
    the equilibrium as the mean field and scale it to obtain a best-fit
    of the station data.
    """
    
    def __init__(self, lock_gamma = None):
        self.lock_gamma = lock_gamma != None
        self.gamma = 1.0 if lock_gamma == None else lock_gamma
        
        # configure diagnostics
        diagnostics().push("mfm_lock_gamma", lock_gamma)
        diagnostics().configure_tag("mfm_res_avg", True, True, True)
        diagnostics().configure_tag("mfm_gamma", True, True, True)
        diagnostics().configure_tag("mfm_mape", True, True, True)
    
    
    def fit_to_data(self, covar, obs, Sigma = None):
        """
        Fit the covariates to the observations using a general least squares
        approach that works for universal kriging (general covariance)
        and for trend surface modelling (diagonal covariance).

        The covariates must already correspond to observation locations.

        Sigma is the covariance matrix of the errors.
        """

        # gather observation data and corresponding model data     
        # FIXME: doing an inverse here instead of a linear solve as I am lazy
        # and I know that matrix is well conditioned
        XtW = np.dot(covar.T, Sigma) if Sigma is not None else covar.T
        print XtW.shape

        XtWX_1 = np.linalg.inv(np.dot(XtW, covar))
        print XtWX_1.shape
    
        # compute the weighted regression
        print obs.shape
        gamma = np.dot(XtWX_1, np.dot(XtW, obs))

        # push diagnostics out
        diagnostics().push("mfm_res_avg", np.mean(obs - gamma * covar))
        diagnostics().push("mfm_gamma", gamma)
        diagnostics().push("mfm_mape", np.mean(np.abs(obs - gamma * covar)))
        
        # if the gamma is not locked, then store the new best fit
        if not self.lock_gamma:
            self.gamma = gamma
            

    def predict_field(self, Covar):
        """
        Return the predicted field given the fit gamma vector.
        The covariates passed must be 3D (lon x lat x covar_id).
        """
        return np.sum(self.gamma[np.newaxis, np.newaxis,:] * Covar, axis = 3)
        
