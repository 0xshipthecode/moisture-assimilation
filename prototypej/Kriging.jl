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

    The kriging results in the matrix K, which contains the kriged observations
    and the matrix V, which contains the kriging variance.
    """
    Nobs = length(obs_data)
    Ncov = size(X,3)

    dsize = size(X)[1:2]
    y = zeros((Nobs,1))
    Xobs = zeros(Nobs, Ncov)

    for (obs,i) in zip(obs_data, 1:Nobs)
    	p = nearest_grid_point(obs)
        y[i] = obs_value(obs)
        Xobs[i,:] = X[p[1], p[2], :]
    end

    # rescale columns of Xobs to have the same norm to reduce
    # the condition number of X'X cheaply
    sc = zeros(Float64, Ncov)
    sc[1] = 1.0
    fm10_norm = sum(Xobs[:,1].^2)^0.5
    for i in 2:Ncov
        sc[i] = fm10_norm / sum(Xobs[:,i].^2)^0.5
        Xobs[:,i] *= sc[i]
    end

    # FIXME: we assume that the measurement variance is the same for
    # all stations at a particular time
    sigma2 = obs_variance(obs_data[1])

    # compute the OLS fit of the covariates to the observations
    XtX = Xobs' * Xobs
    spush("kriging_xtx_cond", cond(XtX))
    beta = XtX \ (Xobs' * y)
    spush("kriging_errors", (Xobs * beta - y)')

    # rescale beta back to the original data (in X)
    beta = beta .* sc
    spush("kriging_beta", beta)

    # compute kriging field and kriging variance and fill out
    # the passed arrays
    for i in 1:dsize[1]
        for j in 1:dsize[2]
            X_ij = squeeze(X[i,j,:], 1)'   # convert covariates at position i,j into a column vector
            K[i,j] = dot(vec(X_ij), vec(beta))
            V[i,j] = sigma2 * (1 + dot(vec(X_ij), vec(XtX \ X_ij)))
        end
    end

end


end