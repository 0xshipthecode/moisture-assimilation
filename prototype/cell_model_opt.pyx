
import numpy as np
cimport numpy as np

cdef extern from "math.h":
    double exp(double x)



cdef class CellMoistureModel:

    cdef np.ndarray Tk                      # nominal fuel delays
    cdef float r0                           # threshold rainfall [mm/h]
    cdef float rk                           # saturation rain intensity [mm/h]
    cdef float Trk                          # time constant for wetting model [s]
    cdef float S                            # saturation intensity [dimensionless]
    
    cdef np.ndarray model_ids, m_ext, m_new, rlag, equi, H, J, change, P2, P
    
    cdef tuple latlon
    
    
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
        else:
            self.Tk = np.array([1, 10, 100]) * 3600.0
            
        self.r0 = 0.05
        self.rk = 8
        self.Trk = 14 * 3600
        self.S = 2.5

        self.model_ids = np.zeros((k,), dtype = np.int)
        self.m_new = np.zeros((k,))
        self.rlag = np.zeros((k,))
        self.equi = np.zeros((k,))
        self.change = np.zeros((k,))
        
        # state covariance matrix
        self.P = np.eye(2*k+3) * 0.02 if P0 is None else P0.copy()
        self.P2 = np.zeros_like(self.P)
        self.H = np.zeros((2*k+3,))
        self.J = np.zeros((2*k+3,2*k+3))

        

    def advance_model(self, np.float64_t Ed,
                            np.float64_t Ew,
                            np.float64_t r,
                            np.float64_t dt,
                            np.ndarray[np.float64_t, ndim=2] mQ = None):
        """
        This model captures the moisture dynamics at each grid point independently.

        Ed - drying equilibrium
        Ew - wetting equilibrium
        r - rain intensity for time unit [mm/h]
        dt - integration step [s]
        """
        
        cdef np.ndarray[np.float64_t, ndim=1] Tk = self.Tk
        cdef int k = Tk.shape[0]                            # number of fuel classes
        cdef float r0 = self.r0
        cdef float S = self.S
        cdef float rk = self.rk
        cdef float Trk = self.Trk

        cdef int i
        
        # first, we break the state vector into components
        cdef np.ndarray[np.float64_t, ndim=1] m_ext = self.m_ext
        
        cdef np.ndarray[np.float64_t, ndim=1] m = m_ext[:k]
        cdef np.ndarray[np.float64_t, ndim=1] dlt_Tk = m_ext[k:2*k]
        cdef np.float64_t dlt_E = m_ext[2*k]
        cdef np.float64_t dlt_S = m_ext[2*k+1]
        cdef np.float64_t dlt_Trk = m_ext[2*k+2]
        
        # add assimilated difference, which is shared across spatial locations
        Ed += dlt_E
        Ew += dlt_E
        
        # where rainfall is above threshold (spatially different), apply
        # saturation model, equi and rlag are specific to fuel type and
        # location
        cdef np.ndarray[np.float64_t, ndim=1] equi = self.equi
        equi[:] = m[:]
        cdef np.ndarray[np.float64_t, ndim=1] rlag = self.rlag
        rlag[:] = 0.0 
        cdef np.ndarray[np.int_t, ndim=1] model_ids = self.model_ids
        model_ids[:] = 0
        
        # equilibrium is equal to the saturation level (assimilated)
        if r > r0:
            equi[:] = S + dlt_S
            model_ids[:] = 3
        
            # rlag is modified by the rainfall intensity (assimilated)
            rlag[:] = 1.0 / (Trk + dlt_Trk) * (1 - exp(- (r - r0) / rk))
            
        else:
            # equilibrium is selected according to current moisture level
            for i in range(k):
                model_ids[i] = 4
                if equi[i] > Ed:
                    model_ids[i] = 1
                    equi[i] = Ed
                elif equi[i] < Ew:
                    model_ids[i] = 2
                    equi[i] = Ew
    
            # the inverted time lag is constant according to fuel category
            rlag[:] = 1.0 / (Tk + dlt_Tk)
        
        # select appropriate integration method according to change for each fuel 
        # and location
        cdef np.ndarray[np.float64_t, ndim=1] change = self.change
        change[:] = rlag[:]
        change *= dt
        cdef np.ndarray[np.float64_t, ndim=1] m_new = self.m_new
        for i in range(k):
            if change[i] < 0.01:
                m_new[i] = m[i] + (equi[i] - m[i]) * (1.0 - exp(-change[i]))
            else:
                m_new[i] = m[i] + (equi[i] - m[i]) * change[i] * (1 - 0.5 * change[i])
        
        # update model state covariance if requested using the old state (the jacobian must be computed as well)
        cdef np.ndarray[np.float64_t, ndim=2] J = self.J
        cdef float dmi_dchng, dmi_dequi 
        if mQ is not None:
            for i in range(k):
            
                if change[i] < 0.01:
                    
                    # partial m_i/partial m_i
                    J[i,i] = exp(-change[i])
                    
                    # precompute partial m_i/partial change
                    dmi_dchng = (equi[i] - m[i]) * exp(-change[i])
                    
                    # precompute partial m_i/partial equi
                    dmi_dequi = (1.0 - exp(-change[i]))
        
                else:
                    
                    # partial dm_i/partial m_i
                    J[i,i] = 1.0 - change[i] * (1 - 0.5 * change[i])
                    
                    # partial m_i/partial change
                    dmi_dchng = (equi[i] - m[i]) * (1.0 - change[i])
                    
                    # partial m_i/partial equi
                    dmi_dequi = change[i] * (1 - 0.5 * change[i])
                                
            
                # branch according to the currently active model
                if r <= r0:
                    
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
                    J[i,2*k+2] = dmi_dchng * dt * (exp(-(r - r0)/self.rk) - 1.0) * (self.Trk + dlt_Trk)**(-2)
        
            
                # delta_Tk for each fuel have no dependencies except previous delta_Tk
                J[k+i,k+i] = 1.0;        
        
            # the equilibrium constants
            J[2*k,2*k] = 1.0 
            J[2*k+1,2*k+1] = 1.0 
            J[2*k+2,2*k+2] = 1.0

            # transformed to run in-place
            np.dot(J, self.P, self.P2)
            np.dot(self.P2, J.T, self.P)
            self.P += mQ

        # update to the new state
        self.m_ext[:3] = m_new


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