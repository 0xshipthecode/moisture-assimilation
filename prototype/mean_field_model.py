
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
        diagnostics().push("mfm_lock_gamma", lock_gamma)
        diagnostics().configure_tag("mfm_res_avg", True, True, True)
        diagnostics().configure_tag("mfm_gamma", True, True, True)
    
    
    def fit_to_data(self, E, obs_list):
        """
        Fit the model to the observations.  Computes a weighted least squares
        estimate of the observed data.
        """
        # gather observation data and corresponding model data        
        grid_pts = [obs.get_nearest_grid_point() for obs in obs_list]
        obsv = np.array([obs.get_value() for obs in obs_list])
        modv = np.array([E[pos] for pos in grid_pts])
        weights = np.array([1.0 / obs.get_station().get_dist_to_grid() for obs in obs_list])
#            weights = np.array([1.0 for obs in obs_list])
    
        # compute the weighted regression
        gamma = np.sum(weights * modv * obsv) / np.sum(weights * modv ** 2)
        diagnostics().push("mfm_res_avg", np.mean(obsv - gamma * modv))
        diagnostics().push("mfm_gamma", gamma)
        
        # if the gamma is not locked, then store the new best fit
        if not self.lock_gamma:
            self.gamma = gamma
            

    def predict_field(self, E):
        """
        Return the predicted field given the fit.
        """
        return self.gamma * E
        
