module DoseResponse

using Distributions
using LinearAlgebra
using Random
using Statistics

"""
    dr_log_rn_fun(dose, rmin, rmax, beta, nec)

Direct translation of `DRLogRNFun` from the Mathematica notebook. Returns the
expected response for a given `dose` using the same piecewise exponential decay
and the same exponential parameterization of `beta`.
"""
function dr_log_rn_fun(dose, rmin, rmax, beta, nec)
    logyhat = rmin + (rmax - rmin) * exp(-exp(beta) * (dose - nec) * (dose > nec))
    return exp(logyhat)
end

"""
    dr_fun_log(dose, alpha, beta, nec)

Literal translation of the piecewise log-response `DRFunLog`.
"""
function dr_fun_log(dose, alpha, beta, nec)
    if dose < nec
        return log(alpha)
    else
        return log(alpha) - exp(beta) * (dose - nec)
    end
end

"""
    dr_fun_log_exp(dose, alpha, beta, nec)

Wrapper matching `DRFunLogExp` that exponentiates the log-response.
"""
dr_fun_log_exp(dose, alpha, beta, nec) = exp(dr_fun_log(dose, alpha, beta, nec))

"""
    generate_individuals(x0, nec, CVx, CVnec, rho, nID; rng = Random.default_rng())

Translation of `GenerateIndividuals`. Draws `nID` correlated `(x0, NEC)` pairs
from a multivariate normal distribution with the provided coefficients of
variation and correlation `rho`.
"""
function generate_individuals(x0, nec, CVx, CVnec, rho, nID; rng = Random.default_rng())
    μ = [x0, nec]
    sigmas = [x0 * CVx, nec * CVnec]
    rho_mat = [1.0 rho; rho 1.0]
    Σ = Diagonal(sigmas) * rho_mat * Diagonal(sigmas)
    dist = MvNormal(μ, Σ)
    samples = permutedims(rand(rng, dist, nID))

    return [
        Dict("x0" => sample[1], "NEC" => sample[2], "ID" => i)
        for (i, sample) in enumerate(eachrow(samples))
    ]
end

"""
    simulate_dose_response(individuals; dose_values, beta, sigma, rng)

Port of `SimulateDoseResponse`. For each individual and dose it draws a
log-normal observation around the deterministic expectation and returns the
per-dose summary statistics (mean and standard deviation of `yhat`).
"""
function simulate_dose_response(individuals; dose_values=0:1:100, beta=-3.5, sigma=0.08, rng=Random.default_rng())
    records = Vector{Dict{String, Float64}}()

    for ind in individuals
        x0 = ind["x0"]
        nec = ind["NEC"]
        id = ind["ID"]
        for dose in dose_values
            logyhat = dr_fun_log(dose, x0, beta, nec)
            yhat = exp(logyhat)
            y = rand(rng, LogNormal(logyhat, sigma))
            push!(records, Dict("ID" => id, "Dose" => float(dose), "logyhat" => logyhat, "yhat" => yhat, "y" => y))
        end
    end

    by_dose = Dict{Float64, Vector{Dict{String, Float64}}}()
    for record in records
        push!(get!(by_dose, record["Dose"], Vector{Dict{String, Float64}}()), record)
    end

    summary = Vector{Dict{String, Float64}}()
    for (dose, vals) in by_dose
        yhat_values = [v["yhat"] for v in vals]
        push!(summary, Dict(
            "Dose" => dose,
            "yhatMean" => mean(yhat_values),
            "yhatSD" => std(yhat_values),
        ))
    end

    sort!(summary, by = d -> d["Dose"])
    return summary
end

end # module
