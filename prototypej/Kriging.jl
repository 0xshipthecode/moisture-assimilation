module Kriging

#
#  The Kriging module provides two types of service:
#
#  * universal kriging with isotropic covariance or correlation
#  * trend surface model kriging
#
#

using Stations
import Stations.nearest_grid_point, Stations.obs_variance, Stations.obs_value

using Storage
import Storage.spush


function universal_kriging(obs, obs_stds, m, m_stds, wrf_data, t)

# Universal kriging not implemented yet

end



function trend_surface_model_kriging(obs_data, X, K, V)
    """
    Trend surface model kriging, which assumes spatially uncorrelated errors.

    WARNING: The variable X is clobbered.

    The kriging results in the matrix K, which contains the kriged observations
    and the matrix V, which contains the kriging variance.
    """
    Nobs = length(obs_data)
    Ncov = size(X,3)

    dsize = size(X)[1:2]
    y = zeros((Nobs,1))
    Xobs = zeros(Nobs, Ncov)
    m_var = zeros(Nobs)

    for (obs,i) in zip(obs_data, 1:Nobs)
    	p = nearest_grid_point(obs)
        y[i] = obs_value(obs)
        Xobs[i,:] = X[p[1], p[2], :]
        m_var[i] = obs_variance(obs)
    end

    # initialize iterative algorithm
    s2_eta_hat_old = -10.0
    s2_eta_hat = 0.0

    i = 0
    while abs(s2_eta_hat_old - s2_eta_hat) > 1e-5
    
        s2_eta_hat = s2_eta_hat_old
        Sigma = diagm(m_var) + s2_eta_hat * eye(Nobs)
        SigInv = inv(Sigma)
        XtSX = Xobs' * SigInv * Xobs
        beta = XtSX \ Xobs' * SigInv * Xobs
        res = y - Xobs * beta
        s2_eta_hat = 1.0 / (Nobs - Ncov) * sum(max(res.^2 - m_var, 0))
        i += 1
    end

    # use the last s2_eta_hat for our estimate
    Sigma = s2_eta_hat * eye(Nobs) + diagm(m_var)
    SigInv = inv(Sigma)

    # compute the OLS fit of the covariates to the observations
    spush("kriging_xtx_cond", cond(XtSX))
    beta = XtSX \ (Xobs' * SigInv * y)
    spush("kriging_errors", (Xobs * beta - y)')

    spush("kriging_beta", beta)

    # compute kriging field and kriging variance and fill out
    # the passed arrays
    for i in 1:dsize[1]
        for j in 1:dsize[2]
            x_ij = squeeze(X[i,j,:], 1)'   # convert covariates at position i,j into a column vector
            K[i,j] = dot(vec(x_ij), vec(beta))
            V[i,j] = s2_eta_hat + dot(vec(x_ij), vec(XtSX \ x_ij))
        end
    end

end


end