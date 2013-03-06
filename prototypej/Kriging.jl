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
    dsize = size(X)[1:2]
    y = zeros((Nobs,1))
    Xobs = zeros((Nobs, size(X,3)))

    for (obs,i) in zip(obs_data, 1:Nobs)
    	p = nearest_grid_point(obs)
        y[i] = obs_value(obs)
        Xobs[i,:] = X[p[1], p[2], :]
    end

    # FIXME: we assume that the measurement variance is the same for
    # all stations at a particular time
    sigma2 = obs_variance(obs_data[1])

    # compute the OLS fit of the covariates to the observations
    XtX = Xobs' * Xobs
    spush("kriging_xtx_cond", cond(XtX))
    beta = XtX \ (Xobs' * y)
    spush("kriging_beta", beta')

    spush("kriging_errors", (Xobs * beta - y)')

    # compute kriging field and kriging variance and fill out
    # the passed arrays
    for i in 1:dsize[1]
        for j in 1:dsize[2]
            X_ij = squeeze(X[i,j,:], 1)'   # convert covariates at position i,j into a column vector
            K[i,j] = (X_ij' * beta)[1,1]
            V[i,j] = sigma2 * (1 + X_ij' * (XtX \ X_ij))[1,1]
        end
    end

end


end