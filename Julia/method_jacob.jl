module MethodJacob

using QuadGK

const INT1_CACHE = Dict{NTuple{10, Float64}, Float64}()
const INT2_CACHE = Dict{NTuple{10, Float64}, Float64}()

pdist(x, X, sigma) = (1 / sqrt(2 * π * sigma)) * exp(-((x - X)^2) / (2 * sigma))
attack(x, a, alpha, τ) = alpha * exp(-((x - a)^2) / (2 * τ^2))
handling(x, h, eta_min, eta_max, nu) = eta_max - (eta_max - eta_min) * exp(-((x - h)^2) / (2 * nu^2))

function _integration_limits(X, sigma)
    if sigma < 1
        return (X - 5, X + 5)
    else
        return (X - 5sigma, X + 5sigma)
    end
end

function int1(R, a, alpha, τ, h, eta_min, eta_max, nu, X, sigma)
    key = (R, a, alpha, τ, h, eta_min, eta_max, nu, X, sigma)
    if haskey(INT1_CACHE, key)
        return INT1_CACHE[key]
    end

    xmin, xmax = _integration_limits(X, sigma)
    integrand(x) = (attack(x, a, alpha, τ) / (1 + attack(x, a, alpha, τ) * handling(x, h, eta_min, eta_max, nu) * R)) * pdist(x, X, sigma)
    value, _ = quadgk(integrand, xmin, xmax; rtol=1e-6, atol=1e-6)
    INT1_CACHE[key] = value
    return value
end

function int2(R, a, alpha, τ, h, eta_min, eta_max, nu, X, sigma)
    key = (R, a, alpha, τ, h, eta_min, eta_max, nu, X, sigma)
    if haskey(INT2_CACHE, key)
        return INT2_CACHE[key]
    end

    xmin, xmax = _integration_limits(X, sigma)
    integrand(x) = begin
        A = attack(x, a, alpha, τ)
        H = handling(x, h, eta_min, eta_max, nu)
        pd = pdist(x, X, sigma)
        (A / (1 + A * H * R)) * pd * (1 - (A * H * R) / (1 + A * H * R))
    end
    value, _ = quadgk(integrand, xmin, xmax; rtol=1e-6, atol=1e-6)
    INT2_CACHE[key] = value
    return value
end

end # module
