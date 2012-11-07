
import numpy as np

from spatial_model_utilities import equilibrium_moisture


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
        self.t = 0.0
        self.model_ids = np.zeros((k,))
        self.m_new = np.zeros((k,))
        self.rlag = np.zeros((k,))
        self.equi = np.zeros((k,))
        
        # state covariance matrix
        self.P = np.eye(2*k+3) * 0.02 if P0 is None else P0.copy()
        self.H = np.zeros((2*k+3,))
        

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
            rlag[:] = 1.0 / (self.Trk + dlt_Trk) * (1 - np.exp(- (r - self.r0) / self.rk))
        else:
            # equilibrium is selected according to current moisture level
            model_ids[:] = 4
            model_ids[equi > Ed] = 1
            equi[equi > Ed] = Ed
            model_ids[equi < Ew] = 2
            equi[equi < Ew] = Ew
    
            # the inverted time lag is constant according to fuel category
            rlag = 1.0 / (Tk + dlt_Tk)
        
        # select appropriate integration method according to change for each fuel 
        # and location
        change = dt * rlag
        m_new = self.m_new
        m_new = np.where(change < 0.01,
                         m + (equi - m) * (1 - np.exp(-change)),
                         m + (equi - m) * change * (1 - 0.5 * change))
        
        # update model state covariance if requested using the old state (the jacobian must be computed as well)
        if mQ is not None:
            J = np.zeros((2*k+3,2*k+3))
            for i in range(k):
            
                if change[i] < 0.01:
                    
                    # partial m_i/partial m_i
                    J[i,i] = np.exp(-change[i])
                    
                    # precompute partial m_i/partial change
                    dmi_dchng = (equi[i] - m[i]) * np.exp(-change[i])
                    
                    # precompute partial m_i/partial equi
                    dmi_dequi = (1.0 - np.exp(-change[i]))
        
                else:
                    
                    # partial dm_i/partial m_i
                    J[i,i] = 1.0 - change[i] * (1 - 0.5 * change[i])
                    
                    # partial m_i/partial change
                    dmi_dchng = (equi[i] - m[i]) * (1.0 - change[i])
                    
                    # partial m_i/partial equi
                    dmi_dequi = change[i] * (1 - 0.5 * change[i])
                                
            
                # branch according to the currently active model
                if r <= self.r0:
                    
                    # drying/wetting model active
        
                    # partial m_i/partial delta_Tk
                    J[i,k+i] = dmi_dchng * (-dt) * (Tk[i] + dlt_Tk[i])**(-2)
        
                    # if drying/wetting model active, jacobian entry w.r.t. equilibrium is nonzero
                    # it is zero if the 'dead zone' model is active
                    if model_ids[i] < 4:
                        J[i, 2*k] = dmi_dequi
        
                else:
        
                    # rain model active
        
                    # partial m_i/partial deltaS
                    J[i,2*k+1] = dmi_dequi
        
                    # partial m_i/partial deltaTkr
                    J[i,2*k+2] = dmi_dchng * dt * (np.exp(-(r - self.r0)/self.rk) - 1.0) * (self.Trk + dlt_Trk)**(-2)
        
            
                # delta_Tk for each fuel have no dependencies except previous delta_Tk
                J[k+i,k+i] = 1.0;        
        
            # the equilibrium constants
            J[2*k,2*k] = 1.0 
            J[2*k+1,2*k+1] = 1.0 
            J[2*k+2,2*k+2] = 1.0

            self.P = np.dot(np.dot(J, self.P), J.T) + mQ

        # update to the new state
        self.m_ext[:3] = m_new
        
        
    def compute_jacobian(self, Ed, Ew, r, dt):
        """
        Compute the jacobian at the current state and return it.
        """
        Tk = self.Tk
        k = Tk.shape[0]                 # number of fuel classes
        
        # first, we break the state vector into components
        m = self.m_ext[:k]
        dlt_Tk = self.m_ext[k:k+k]
        dlt_E = self.m_ext[2*k]
        dlt_S = self.m_ext[2*k+1]
        dlt_Trk = self.m_ext[2*k+2]
        
        # add assimilated difference, which is shared across spatial locations
        Ed = Ed + dlt_E
        Ew = Ew + dlt_E
    
        # where rainfall is above threshold (spatially different), apply
        # saturation model, equi and rlag are specific to fuel type and
        # location
        equi = m.copy()         # copy over current equilibrium levels
        rlag = np.zeros((k,))
        model_ids = np.zeros((k,))
        
        # equilibrium is equal to the saturation level (assimilated)
        if r > self.r0:
            equi[:] = self.S + dlt_S
            self.model_ids[:] = 3
        
            # rlag is modified by the rainfall intensity (assimilated)
            rlag[:] = 1.0 / (self.Trk + dlt_Trk) * (1 - np.exp(- (r - self.r0) / self.rk))
        else:
    
            # equilibrium is selected according to current moisture level
            self.model_ids[:] = 4
            self.model_ids[equi > Ed] = 1
            equi[equi > Ed] = Ed
            self.model_ids[equi < Ew] = 2
            equi[equi < Ew] = Ew
    
            # the inverted time lag is constant according to fuel category
            rlag = 1.0 / (Tk + dlt_Tk)

        # select appropriate integration method according to change for each fuel 
        # and location
        change = dt * rlag
        J = np.zeros((2*k+3,2*k+3))
        for i in range(k):
        
            if change[i] < 0.01:
                
                # partial m_i/partial m_i
                J[i,i] = np.exp(-change[i])
                
                # precompute partial m_i/partial change
                dmi_dchng = (equi[i] - m[i]) * np.exp(-change[i])
                
                # precompute partial m_i/partial equi
                dmi_dequi = (1.0 - np.exp(-change[i]))
    
            else:
                
                # partial dm_i/partial m_i
                J[i,i] = 1.0 - change[i] * (1 - 0.5 * change[i])
                
                # partial m_i/partial change
                dmi_dchng = (equi[i] - m[i]) * (1.0 - change[i])
                
                # partial m_i/partial equi
                dmi_dequi = change[i] * (1 - 0.5 * change[i])
                            
        
            # branch according to the currently active model
            if r <= self.r0:
                
                # drying/wetting model active
    
                # partial m_i/partial delta_Tk
                J[i,k+i] = dmi_dchng * (-dt) * (Tk[i] + dlt_Tk[i])**(-2)
    
                # if drying/wetting model active, jacobian entry w.r.t. equilibrium is nonzero
                # it is zero if the 'dead zone' model is active
                if model_ids[i] < 4:
                    J[i, 2*k] = dmi_dequi
    
            else:
    
                # rain model active
    
                # partial m_i/partial deltaS
                J[i,2*k+1] = dmi_dequi
    
                # partial m_i/partial deltaTkr
                J[i,2*k+2] = dmi_dchng * dt * (np.exp(-(r - self.r0)/self.rk) - 1.0) * (self.Trk + dlt_Trk)**(-2)
    
        
            # delta_Tk for each fuel have no dependencies except previous delta_Tk
            J[k+i,k+i] = 1.0;        
    
        # the equilibrium constants
        J[2*k,2*k] = 1.0 
        J[2*k+1,2*k+1] = 1.0 
        J[2*k+2,2*k+2] = 1.0
        
        return J 


    def get_state(self):
        """
        Return a copy of the current state.
        """
        return self.m_ext.copy()
    
    
    def kalman_update(self, obs_val, obs_var, fuel_type):
        """
        Updates the state using the observation at the grid point.
        
          obs_val - the value of the observed fuel moisture
          obs_var - the variance of the observed fuel moisture
          fuel_type - the fuel type [0..k-1], where k is the number of fuels run in the model
        """
        
        P, H = self.P, self.H
        H[:] = 0.0
        H[fuel_type] = 1.0
        
        # we assume we only observe one fuel at a time (simplifes inversion)
        I = np.dot(np.dot(H, P), H.T) + obs_var
        K = np.dot(P, H.T) / I
        
        # update state and state covariance
        self.m_ext += K * (obs_val - self.m_ext[fuel_type])
        P -= np.dot(np.dot(K, H), P) 

        return K[fuel_type]