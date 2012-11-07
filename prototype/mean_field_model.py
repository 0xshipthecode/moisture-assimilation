
import numpy as np


class MeanFieldModel:
    """
    The model of the mean field for simple kriging.  The idea is to take
    the equilibrium as the mean field and scale it to obtain a best-fit
    of the station data.
    """
    
    def __init__(self):
        self.gamma = 1.0
    
    
    def fit_to_data(self, E, obs_list):
        """
        Fit the model to the observations.  Computes a weighted least squares
        estimate of the observed data.  The weight is given by the inverse
        distance from the nearest grid point to the observation station.
        """
        # gather observation data and corresponding model data        
        grid_pts = [obs.get_nearest_grid_point() for obs in obs_list]
        obsv = np.array([obs.get_value() for obs in obs_list])
        modv = np.array([E[pos] for pos in grid_pts])
        weights = np.array([1.0 / obs.get_station().get_dist_to_grid() for obs in obs_list])
        
        # compute the weighted regression
        self.gamma = np.sum(weights * modv * obsv) / np.sum(weights * modv ** 2)


    def predict_field(self, E):
        """
        Return the predicted field given the fit.
        """
        return self.gamma * E
        