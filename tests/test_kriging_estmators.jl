

function trend_surface_model_kriging(X, Z, Sigma2_eps)

    # now run the iterative kriging estimator
    sigma2_eta_hat = 0
    sigma2_eta_hat_old = -1
    beta_hat = []
    res = []
    i = 1
    while abs(sigma2_eta_hat - sigma2_eta_hat_old) > 1e-4
        
        sigma2_eta_hat_old = sigma2_eta_hat
#        println("Iteration $i:")
        
        # construct Sigma matrix based on our best estimate
        Sigma = Sigma2_eps + sigma2_eta_hat * eye(Nobs)

        # estimate beta first
        SigInv = diagm(1.0 ./ diag(Sigma))
        beta_hat = inv(X' * SigInv * X) * X' * SigInv * Z

        # compute the "kriging residuals"
        res = Z - X * beta_hat
#        println("Residual mean $(mean(res))")

        # then estimate sigma_eta
        sigma2_eta_hat = 1.0 / (size(X,1) - size(X,2)) * sum(max(res.^2 - diag(Sigma2_eps), 0))

#        println("  beta_hat $beta_hat")
#        println("  sigma2_eta_hat $sigma2_eta_hat")

        i += 1
    end

    return sigma2_eta_hat, res

end


function trend_surface_model_kriging_bootstrap(X, Z, Sigma2_eps)

    Nobs = size(X,1)
    Niters = 500
    s2_etas = zeros(Niters)
    res_sum = zeros(Nobs, 1)

    s2_eta_hat_orig, res_i = trend_surface_model_kriging(X, Z, Sigma2_eps)

    for i=1:Niters

        ndx = rand(1:Nobs,Nobs)
        
        Xi = X[ndx, :]
        Zi = Z[ndx,:]
        Sigma2_epsi = Sigma2_eps[ndx,ndx]

        s2_eta_hat, res_i = trend_surface_model_kriging(Xi, Zi, Sigma2_epsi)
        s2_etas[i] = s2_eta_hat
        res_sum += res_i

    end

#    s2_eta_hat_mean = 1.0 / (size(X,1) - size(X,2)) * sum(1.0/Niters * res_sum - diag(Sigma2_eps))
    s2_eta_hat_mean = sum(s2_etas) / Niters

    return s2_eta_hat_orig * 2 - s2_eta_hat_mean, s2_etas

end


# first construct an actual model
Nobs = 20
beta = [ 1.1, -0.2, 0.3 ]''
Sigma2_eps = 0.02^2 * eye(Nobs)
sigma2_eta = 0.04^2

s2 = Float64[]
s2_etas = nothing

for i = 1:200

    X = randn(Nobs,3)
    Z_exact = X * beta
    Z = X * beta + sqrt(Sigma2_eps) * randn((Nobs,1)) + sqrt(sigma2_eta) * randn((Nobs,1))

#    sigma2_eta_hat, res_i = trend_surface_model_kriging_bootstrap(X, Z, Sigma2_eps)
    sigma2_eta_hat, s2_etas = trend_surface_model_kriging_bootstrap(X, Z, Sigma2_eps)
    push!(s2, sigma2_eta_hat)

end

println(sigma2_eta)
#println(sigma2_eta_hat)
for s in s2
    println(s)
end

    