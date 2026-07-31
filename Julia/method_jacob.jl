using QuadGK
using Roots
using LinearAlgebra

# Probability distribution
pdist(y, Y, sigma) = (1 / sqrt(2π * sigma)) * exp(-(y - Y)^2 / (2sigma))

# Attack rate
attack(y, a, alpha, tau) = alpha * exp(-(y - a)^2 / (2tau^2))

# Handling rate
handling(y, h, eta_min, eta_max, nu) = eta_max .- (eta_max .- eta_min) * exp(-(y - h)^2 / (2nu^2))

# Parameter values (translated from Mathematica `params_gen`)
params_gen = (
    r       = 0.3,
    alpha   = 2,
    eta_max = 2,
    eta_min = 1,
    epsilon = 0.5,
    tau     = 1,
    nu      = 1,
    k       = 1,
    beta    = 0.1,
    Y       = 3,
    a       = 2.5,
    h       = 2.5,
    sigma   = sqrt(3)
)

# Integral 1 (equivalent to Int1)
function int_1(R, a, alpha, tau, h, eta_min, eta_max, nu, Y, sigma)
    f(y) = (attack(y, a, alpha, tau) /
           (1 + attack(y, a, alpha, tau) * handling(y, h, eta_min, eta_max, nu) * R)) *
           pdist(y, Y, sigma)
    val, _ = quadgk(f, -Inf, Inf) # This evaluates the integral of the function f(y) from -Inf to +Inf
    return val
end

# Integral 2 (equivalent to Int2)
function int_2(R, a, alpha, tau, h, eta_min, eta_max, nu, Y, sigma)
    f(y) = begin
        atk = attack(y, a, alpha, tau)
        hnd = handling(y, h, eta_min, eta_max, nu)
        denom = 1 + atk * hnd * R
        pd = pdist(y, Y, sigma)
        term = 1 - (atk * hnd * R) / denom
        (atk / denom) * pd * term # This is the value of f(y)
    end
    val, _ = quadgk(f, -Inf, Inf)
    return val
end


# Function whose root we solve for (FNum)
function f_num(R, a, alpha, tau, h, eta_min, eta_max, nu, Y, sigma, beta, epsilon)
    R * int_1(R, a, alpha, tau, h, eta_min, eta_max, nu, Y, sigma) - beta / epsilon
end

# # Solve for R (equivalent to FindRoot)
# R_sol = find_zero(
#     R -> f_num(R,
#                params_gen.a, params_gen.alpha, params_gen.tau,
#                params_gen.h, params_gen.eta_min, params_gen.eta_max, params_gen.nu,
#                params_gen.Y, params_gen.sigma,
#                params_gen.beta, params_gen.epsilon),
#     0.5
# )

# # Compute C_sol (equivalent to C_sol)
# C_sol = (params_gen.r * R_sol * (1 - R_sol / params_gen.k)) /
#        (R_sol * int_1(R_sol,
#                      params_gen.a, params_gen.alpha, params_gen.tau,
#                      params_gen.h, params_gen.eta_min, params_gen.eta_max, params_gen.nu,
#                      params_gen.Y, params_gen.sigma)
#        )

# # println("R_sol = ", R_sol)
# # println("C_sol = ", C_sol)

# j11 = params_gen.r * (1 - (2 * R_sol) / params_gen.k) -
#        C_sol * int_2(R_sol, params_gen.a, params_gen.alpha, params_gen.tau,
#                     params_gen.h, params_gen.eta_min, params_gen.eta_max, params_gen.nu,
#                     params_gen.Y, params_gen.sigma)

# j12 = -R_sol * int_1(R_sol, params_gen.a, params_gen.alpha, params_gen.tau,
#                     params_gen.h, params_gen.eta_min, params_gen.eta_max, params_gen.nu,
#                     params_gen.Y, params_gen.sigma)

# j21 = params_gen.epsilon * C_sol * int_2(R_sol, params_gen.a, params_gen.alpha, params_gen.tau,
#                                         params_gen.h, params_gen.eta_min, params_gen.eta_max, params_gen.nu,
#                                         params_gen.Y, params_gen.sigma)

# j22 = params_gen.epsilon * R_sol * int_1(R_sol, params_gen.a, params_gen.alpha, params_gen.tau,
#                                         params_gen.h, params_gen.eta_min, params_gen.eta_max, params_gen.nu,
#                                         params_gen.Y, params_gen.sigma) - params_gen.beta

# jacobian_at_eq = [
#     j11  j12;
#     j21  j22
# ]

# # println("Jacobian at equilibrium:")
# # println(jacobian_at_eq)

# det_val = det(jacobian_at_eq)
# tr_val = tr(jacobian_at_eq)

# # println("det = ", det_val)
# # println("tr  = ", tr_val)
# # println("Delta = ", tr_val^2 - 4*det_val)

# Function finding the numerical equilibrium (EqSol):
function eq_sol(d_val, sigma_val)
    params_spec = merge(
        params_gen,
        (
            a     = params_gen.Y - d_val,
            h     = params_gen.Y - d_val,
            sigma = sigma_val
        )
    )
    f_num_val(R) = f_num(
        R,
        params_spec.a,
        params_spec.alpha,
        params_spec.tau,
        params_spec.h,
        params_spec.eta_min,
        params_spec.eta_max,
        params_spec.nu,
        params_spec.Y,
        params_spec.sigma,
        params_spec.beta,
        params_spec.epsilon
    )
    R_sol = find_zero(f_num_val, 0.5)
    C_sol = (
        params_spec.r * R_sol * (1 - R_sol / params_spec.k)) / (
            R_sol * int_1(
                R_sol,
                params_spec.a, params_spec.alpha, params_spec.tau,
                params_spec.h, params_spec.eta_min, params_spec.eta_max, params_spec.nu,
                params_spec.Y, params_spec.sigma
            )
        )
    return R_sol, C_sol
end

# test_sol = eq_sol(0.5,sqrt(3))
# println("test_sol = ", test_sol)

# Function computing the det and tr from the equilibrium (DetTr)
function det_tr(d_val, sigma_val, sol)
    R_sol, C_sol = sol
    params_spec = merge(
        params_gen,
        (
            a     = params_gen.Y - d_val,
            h     = params_gen.Y - d_val,
            sigma = sigma_val
        )
    )
    j11 = params_spec.r * (1 - (2 * R_sol) / params_spec.k) -
        C_sol * int_2(
            R_sol, params_spec.a, params_spec.alpha, params_spec.tau,
            params_spec.h, params_spec.eta_min, params_spec.eta_max, params_spec.nu,
            params_spec.Y, params_spec.sigma
        )
    j12 = -R_sol * int_1(
        R_sol, params_spec.a, params_spec.alpha, params_spec.tau,
        params_spec.h, params_spec.eta_min, params_spec.eta_max, params_spec.nu,
        params_spec.Y, params_spec.sigma
    )
    j21 = params_spec.epsilon * C_sol * int_2(
        R_sol, params_spec.a, params_spec.alpha, params_spec.tau,
        params_spec.h, params_spec.eta_min, params_spec.eta_max, params_spec.nu,
        params_spec.Y, params_spec.sigma
    )
    j22 = params_spec.epsilon * R_sol * int_1(
        R_sol, params_spec.a, params_spec.alpha, params_spec.tau,
        params_spec.h, params_spec.eta_min, params_spec.eta_max, params_spec.nu,
        params_spec.Y, params_spec.sigma
    ) - params_spec.beta
    jacobian_at_eq = [
        j11  j12;
        j21  j22
    ]
    det_val = det(jacobian_at_eq)
    tr_val = tr(jacobian_at_eq)
    return det_val, tr_val
end

# test_jacob = det_tr(0.5, sqrt(3), test_sol)
# println("test_jacob = ", test_jacob)

# Function classifying the stability class (Classifying)
function classifying(d_val, sigma_val)
    sol = eq_sol(d_val, sigma_val)  # returns (R, c)
    R_sol, c_sol = sol       # unpack tuple

    if R_sol > 0 && c_sol > 0
        # True: coexistence
        coeffs = det_tr(d_val, sigma_val, sol)  # returns (det, tr)
        det_val, tr_val = coeffs         # unpack

        if det_val > 0 && tr_val < 0
            # True: no limit cycles
            if tr_val^2 - 4 * det_val >= 0
                # True: nonoscillatory
                return 1
            else
                # False: damped oscillations
                return 2
            end
        else
            # False: limit cycles
            return 3
        end
    else
        # False: no coexistence
        return 0
    end
end

# test_class = classifying(0.5, sqrt(3))
# println("stability class: ", test_class)