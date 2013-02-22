
import numpy as np



def construct_covariance_matrix_lonlat(lonlat1, sig1, n, s, lonlat2, sig2):
    """
    Construct a covariance matrix of size len(sig_obs) x len(sig_mod) that 
    is np.diag(sig1) * corr(lonlat) * np.diag(sig2).  In this case n, s
    specify the nugget/slope of the correlation which is converted to covariance.

    sig2/lonlat2 is replaced by sig1/lonlat1 if not given to construct a square matrix.

    If neither sig1 is not given, then it is assumed that the n and s specify covariance
    directly.
    """

    # first construct the correlation matrix
    C = construct_correlation_matrix_lonlat(lonlat1, n, s, lonlat2)

    # now multiply by sig1 from the left
    if sig1 is not None:
        C *= sig1[:, np.newaxis]

        # multiply by sig2 from the right if available, else use sig1
        if sig2 is not None:
            C *= sig2[np.newaxis, :]
        else:
            C *= sig1[np.newaxis, :]

    return C


def construct_correlation_matrix_lonlat(lonlat1, n, s, lonlat2 = None):
    """
    Construct a distance-based correlation matrix between residuals at given longitudes
    and lattitudes.  The correlation structure is isotropic with a nugget term ''n''
    and a slope ''s''.
    """
    N1 = len(lonlat1)
    N2 = N1 if lonlat2 is None else len(lonlat2)
    lonlat2 = lonlat2 if lonlat2 is None

    C = np.zeros((N1, N2))
    
    # compute distances in km between locations
    for (lon1,lat1), i1 in zip(lonlat1, range(N1)):
        for (lon2,lat2), i2 in zip(lonlat2, range(N2)):
            if i1 != i2:
                C[i1,i2] = max(0.0, n - s * great_circle_distance(lon1, lat1, lon2, lat2))
            else:
                C[i1,i2] = 0.0
                
    # rectangular correlation matrix
    return C


def construct_diagonal_covariance(sigmas):
    """
    Construct a diagonal covariance matrix that only contains measurement variances
    for individual stations.
    """
    return np.diag(sigmas)



