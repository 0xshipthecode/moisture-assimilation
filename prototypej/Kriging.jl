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

    # quick pre-conditioning hack
    # rescale all X[:,:,i] to have norm of X[:,:,1]
    n1 = sum(X[:,:,1].^2)^0.5
    for i in 2:Ncov
        X[:,:,i] *= n1 / sum(X[:,:,i].^2)^0.5
    end

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
    subzeros = 0
    while abs(s2_eta_hat_old - s2_eta_hat) > 1e-4
    
        s2_eta_hat_old = s2_eta_hat
        Sigma = diagm(m_var) + s2_eta_hat * eye(Nobs)
        SigInv = inv(Sigma)
        XtSX = Xobs' * SigInv * Xobs
        beta = XtSX \ Xobs' * SigInv * y
        res = y - Xobs * beta
        s2_eta_hat = 1.0 / (Nobs - Ncov) * sum(max(res.^2 - m_var, 0))
        subzeros = sum(res.^2 - m_var .< 0)
        i += 1
        println("Iter: $i  old $s2_eta_hat_old  new $s2_eta_hat")
    end

    # construct new covariance with last s2_eta_hat
    Sigma = s2_eta_hat * eye(Nobs) + diagm(m_var)
    SigInv = inv(Sigma)

    # estimate the beta trend parameters
    XtSX = Xobs' * SigInv * Xobs
    beta = XtSX \ (Xobs' * SigInv * y)

    # compute the OLS fit of the covariates to the observations
    spush("kriging_xtx_cond", cond(XtSX))
    spush("kriging_errors", (Xobs * beta - y)')
    spush("kriging_beta", beta)
    spush("kriging_sigma2_eta", s2_eta_hat)
    spush("kriging_iters", i)
    spush("kriging_subzero_s2_estimates", subzeros)

    # compute kriging field and kriging variance 
    for i in 1:dsize[1]
        for j in 1:dsize[2]
            x_ij = squeeze(X[i,j,:], 1)'   # convert covariates at position i,j into a column vector
            K[i,j] = dot(vec(x_ij), vec(beta))
            V[i,j] = s2_eta_hat + dot(vec(x_ij), vec(XtSX \ x_ij))
        end
    end

end


end