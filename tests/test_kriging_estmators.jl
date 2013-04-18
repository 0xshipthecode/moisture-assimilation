

function numerical_solve_bisect(e2, eps2, k)

    N = size(e2,1)
    tgt = N - k
    s2_eta_left = 0.0
    s2_eta_right = 0.1

    val_left = sum(e2 ./ (eps2 + s2_eta_left))
    if val_left < tgt
        return -1.0
    end

    val_right = sum(e2 ./ (eps2 + s2_eta_right))
    while val_right > tgt
        s2_eta_right *= 2.0
        val_right = sum(e2 ./ (eps2 + s2_eta_right))
    end

    # Newtons method implementation (initialized with s2_eta = 0)
    while abs(val_left - val_right) / (0.5 * (val_right + val_left)) > 1e-3
        
        # compute new value at center of eta interval
        s2_eta = 0.5 * (s2_eta_left + s2_eta_right)
        val = sum(e2 ./ (eps2 + s2_eta))

        if val > tgt
            val_left, s2_eta_left = val, s2_eta
        else
            val_right, s2_eta_right = val, s2_eta
        end

    end

    return 0.5 * (s2_eta_left + s2_eta_right)

end


function trend_surface_model_kriging(X, Z, sigma2_eps, meth)

#    println("***** TSM *****")

    # now run the iterative kriging estimator
    sigma2_eta_hat = 0
    sigma2_eta_hat_old = -1
    beta_hat = []
    res = []
    i = 1
    while abs(sigma2_eta_hat - sigma2_eta_hat_old) / sigma2_eta_hat > 1e-3
        
        sigma2_eta_hat_old = sigma2_eta_hat
        
        # construct Sigma matrix based on our best estimate
        Sigma = (sigma2_eps + sigma2_eta_hat) * eye(Nobs)

        # estimate beta first
        XSX = (X' * (Sigma \ X))
        beta_hat = XSX \ X' * (Sigma \ Z)

        # compute the "kriging residuals"
        res = Z - X * beta_hat

        # then estimate sigma_eta
        s2_eps = sigma2_eps * ones(Float64, size(res))
        if meth == 1
            sigma2_eta_hat = numerical_solve_bisect(res.^2, s2_eps, size(X,2))
        elseif meth == 2
            sigma2_eta_hat = 0.0
            for j=1:size(X,1)
                sigma2_eta_hat += res[j]^2 - s2_eps[j] + dot(vec(X[j,:]), vec((XSX \ X[j,:]')))
            end
            sigma2_eta_hat /= size(X,1)
        elseif meth == 3
            sigma2_eta_hat = sum(res.^2) / (size(X,1) - size(X,2)) - s2_eps[1]
        else
            error("Error in method specification")
        end
            

#        println("iter $i:  beta_hat $beta_hat")
#        println("iter $i:  sigma2_eta_hat $sigma2_eta_hat")

        i += 1
    end

    return sigma2_eta_hat, beta_hat, res

end



# first construct an actual model
Nobs = int(ARGS[1])
meth = int(ARGS[2])

beta = [ 1.1, -0.3 ]''
sigma_eps = 0.02
sigma2_eps = sigma_eps^2
sigma_eta = 0.04
sigma2_eta = sigma_eta^2

p = [beta, sigma2_eta]'
println("$(p[1]), $(p[2]), $(p[3])")

s2 = Array{Float64}[]
s2_etas = nothing
errs = 0

i = 1
while i < 5000

    X = randn(Nobs,2)
    Z_exact = X * beta
    Z = X * beta + sigma_eps * randn((Nobs,1)) + sigma_eta * randn((Nobs,1))

    sigma2_eta_hat, beta_hat, res = trend_surface_model_kriging(X, Z, sigma2_eps, meth)

    if sigma2_eta_hat > 0
        p = [beta_hat, sigma2_eta_hat]'
        println("$(p[1]), $(p[2]), $(p[3])")
        i += 1
    else
        errs += 1
    end

end

println("$errs, 0, 0")
