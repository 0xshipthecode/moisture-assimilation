
import numpy as np
import math


class CellMoistureModel:

    Tk = np.array([1, 10, 100]) * 3600.0    # nominal fuel delays
    r0 = 0.05                               # threshold rainfall [mm/h]
    rk = 8                                  # saturation rain intensity [mm/h]
    Trk = 14 * 3600                         # time constant for wetting model [s]
    S = 2.5                                 # saturation intensity [dimensionless]
    
    
    def __init__(self, latlon, k, m0 = None, Tk = None, P0 = None):
        """
        Initialize the model with given position and moisture levels. 
        """
        self.latlon = latlon
        self.m_ext = np.zeros((2*k+3,))
        if m0 is not None:
            self.m_ext[:k] = m0
        if Tk is not None:
            self.Tk = Tk

        self.model_ids = np.zeros((k,))
        self.m_new = np.zeros((k,))
        self.rlag = np.zeros((k,))
        self.equi = np.zeros((k,))
        
        # state covariance matrix
        self.P = np.eye(2*k+3) * 0.02 if P0 is None else P0.copy()
        self.P2 = np.zeros_like(self.P)
        self.H = np.zeros((k, 2*k+3))
        self.J = np.zeros((2*k+3, 2*k+3))
        

    def advance_model(self, Ed, Ew, r, dt, mQ = None):
        """
        This model captures the moisture dynamics at each grid point independently.

        Ed - drying equilibrium
        Ew - wetting equilibrium
        r - rain intensity for time unit [mm/h]
        dt - integration step [s]
        """
        Tk = self.Tk
        k = Tk.shape[0]                 # number of fuel classes
        
        # first, we break the state vector into components
        m = self.m_ext[:k]
        dlt_Tk = self.m_ext[k:2*k]
        dlt_E = self.m_ext[2*k]
        dlt_S = self.m_ext[2*k+1]
        dlt_Trk = self.m_ext[2*k+2]
        
        # add assimilated difference, which is shared across spatial locations
        Ed = Ed + dlt_E
        Ew = Ew + dlt_E
        
        # where rainfall is above threshold (spatially different), apply
        # saturation model, equi and rlag are specific to fuel type and
        # location
        equi = self.equi
        equi[:] = m[:]
        rlag = self.rlag
        rlag[:] = 0.0 
        model_ids = self.model_ids
        model_ids[:] = 0
        
        # equilibrium is equal to the saturation level (assimilated)
        if r > self.r0:
            equi[:] = self.S + dlt_S
            model_ids[:] = 3
        
            # rlag is modified by the rainfall intensity (assimilated)
            rlag[:] = 1.0 / (self.Trk + dlt_Trk) * (1.0 - math.exp(- (r - self.r0) / self.rk))
            
        else:
            # equilibrium is selected according to current moisture level
            model_ids[:] = 4
            model_ids[equi > Ed] = 1
            equi[equi > Ed] = Ed
            model_ids[equi < Ew] = 2
            equi[equi < Ew] = Ew
    
            # the inverted time lag is constant according to fuel category
            rlag[:] = 1.0 / (Tk + dlt_Tk)
        
        # select appropriate integration method according to change for each fuel 
        # and location
        change = dt * rlag
        m_new = self.m_new
        for i in range(k):
            if change[i] < 0.01:
                m_new[i] = m[i] + (equi[i] - m[i]) * (1.0 - math.exp(-change[i]))
            else:
                m_new[i] = m[i] + (equi[i] - m[i]) * change[i] * (1 - 0.5 * change[i])
        
        # update model state covariance if requested using the old state (the jacobian must be computed as well)
        if mQ is not None:
            J = self.J
            J[:] = 0.0
            
            for i in range(k):
            
                if change[i] < 0.01:
                    
                    # partial m_i/partial m_i
                    J[i,i] = math.exp(-change[i]) if model_ids[i] != 4 else 1.0
                    
                    # precompute partial m_i/partial change
                    dmi_dchng = (equi[i] - m[i]) * math.exp(-change[i])
                    
                    # precompute partial m_i/partial equi
                    dmi_dequi = (1.0 - math.exp(-change[i]))
        
                else:
                    
                    # partial dm_i/partial m_i
                    J[i,i] = 1.0 - change[i] * (1 - 0.5 * change[i]) if model_ids[i] != 4 else 1.0
                    
                    # partial m_i/partial change
                    dmi_dchng = (equi[i] - m[i]) * (1.0 - change[i])
                    
                    # partial m_i/partial equi
                    dmi_dequi = change[i] * (1 - 0.5 * change[i])
                                
            
                # branch according to the currently active model
                if r <= self.r0 and model_ids[i] != 4:
                    
                    # partial m_i/delta E
                    J[i, 2*k] = dmi_dequi
        
                    # partial m_i/partial delta_Tk
                    J[i,k+i] = dmi_dchng * (-dt) * (Tk[i] + dlt_Tk[i])**(-2)

                else:
        
                    # rain model active
        
                    # partial m_i/partial deltaS
                    J[i,2*k+1] = dmi_dequi
        
                    # partial m_i/partial deltaTkr
                    J[i,2*k+2] = dmi_dchng * dt * (math.exp(-(r - self.r0)/self.rk) - 1.0) * (self.Trk + dlt_Trk)**(-2)
        
            
                # delta_Tk for each fuel have no dependencies except previous delta_Tk
                J[k+i,k+i] = 1.0;        
        
            # the equilibrium constants
            J[2*k,2*k] = 1.0 
            J[2*k+1,2*k+1] = 1.0 
            J[2*k+2,2*k+2] = 1.0

            # transformed to run in-place with one pre-allocated temporary
            np.dot(J, self.P, self.P2)
            np.dot(self.P2, J.T, self.P)
            self.P += mQ

        # update to the new state
        self.m_ext[:k] = m_new
        
        
    def get_state(self):
        """
        Return the current state. READ-ONLY under normal circumstances.
        """
        return self.m_ext
    
    
    def get_state_covar(self):
        """
        Return the state covariance. READ-ONLY under normal circumstances.
        """
        return self.P
    
    
    def kalman_update(self, O, V, fuel_types):
        """
        Updates the state using the observation at the grid point.
        
          O - the 1D vector of observations 
          V - the (diagonal) variance matrix of the measurements
          fuel_types - the fuel types for which the observations exist (used to construct observation vector)
          
        """
        P, H = self.P, self.H
        H[:] = 0.0
        
        # construct an observation operator tailored to observed fuel types 
        for i in fuel_types:
            H[i,i] = 1.0
        Ho = H[fuel_types,:]
        
        # use general observation model
        I = np.dot(np.dot(Ho, P), Ho.T) + V
        K = np.dot(np.dot(P, Ho.T), np.linalg.inv(I))
        
        # update state and state covariance
        self.m_ext += np.dot(K, O - self.m_ext[fuel_types])
        P -= np.dot(np.dot(K, Ho), P) 

        return K
