module DoseResponse

using Distributions
using LinearAlgebra
using Random
using Statistics
using Plots
using StatsPlots
gr()

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
    mu = [x0, nec]
    sigmas = [x0 * CVx, nec * CVnec]
    rho_mat = [1.0 rho; rho 1.0]
    sigma_mat = Diagonal(sigmas) * rho_mat * Diagonal(sigmas)
    dist = MvNormal(mu, sigma_mat)
    samples = permutedims(rand(rng, dist, nID))

    return [
        Dict("x0" => sample[1], "NEC" => sample[2], "ID" => i)
        for (i, sample) in enumerate(eachrow(samples))
    ]
end

x_0 = 100
beta = -3.5
dose = range(0, 100)
NEC = 25
n_id = 10_000
CV_x = 0.2
CV_NEC = 0.2
sigma = 0.08
rho = 0.9

individuals_pos = generate_individuals(x_0, NEC, CV_x, CV_NEC, rho, n_id)
individuals_neg = generate_individuals(x_0, NEC, CV_x, CV_NEC, -rho, n_id)
individuals_null = generate_individuals(x_0, NEC, CV_x, CV_NEC, 0, n_id)

col_pos = :red
col_neg = :blue
col_null = :green

function plot_distrib(individualsPos, individualsNeg, individualsNull,
                     colPos, colNeg, colNull)

    # Extract {x0, NEC} pairs for scatter plot
    dataPos  = [(ind["x0"], ind["NEC"])  for ind in individualsPos]
    dataNeg  = [(ind["x0"], ind["NEC"])  for ind in individualsNeg]
    dataNull = [(ind["x0"], ind["NEC"]) for ind in individualsNull]

    # ---- Scatter plot (ListPlot analogue) ----
    plt_scatter = scatter(
        first.(dataPos),  last.(dataPos),
        color = colPos, ms = 0.8,
        markerstrokecolor = :transparent,
        markerstrokewidth = 0,
        label = "Positive cov",
        xlabel = "x₀", ylabel = "NEC",
        ylim = (0, 45),
        aspect_ratio = 1, size = (600, 600),
        framestyle=:box, legend=:topright
    )

    scatter!(
        first.(dataNeg), last.(dataNeg),
        color = colNeg, ms = 0.8,
        markerstrokecolor = :transparent,
        markerstrokewidth = 0,
        label = "Negative cov"
    )

    scatter!(
        first.(dataNull), last.(dataNull),
        color = colNull, ms = 0.8,
        markerstrokecolor = :transparent,
        markerstrokewidth = 0,
        label = "No cov"
    )

    # ---- SmoothHistogram analogue using density() ----

    # Extract components for KDE
    pos_x0  = [ind["x0"]  for ind in individualsPos]
    pos_NEC = [ind["NEC"] for ind in individualsPos]

    neg_x0  = [ind["x0"]  for ind in individualsNeg]
    neg_NEC = [ind["NEC"] for ind in individualsNeg]

    null_x0 = [ind["x0"]  for ind in individualsNull]
    null_NEC = [ind["NEC"] for ind in individualsNull]

    plt_density = plot(
        xlabel="Value", ylabel="Density",
        framestyle=:box, legend=:topright,
        title="Smoothed Distributions"
    )

    # Positive
    density!(plt_density, pos_x0,  color=colPos,  lw=2, label="Positive cov (x₀)")
    density!(plt_density, pos_NEC, color=colPos,  lw=2, linestyle=:dash, label="Positive cov (NEC)")

    # Negative
    density!(plt_density, neg_x0,  color=colNeg,  lw=2, label="Negative cov (x₀)")
    density!(plt_density, neg_NEC, color=colNeg,  lw=2, linestyle=:dash, label="Negative cov (NEC)")

    # Null
    density!(plt_density, null_x0,  color=colNull, lw=2, label="Null cov (x₀)")
    density!(plt_density, null_NEC, color=colNull, lw=2, linestyle=:dash, label="Null cov (NEC)")

    return plt_scatter, plt_density
end

scatter_plot, hist_plot = plot_distrib(individuals_pos, individuals_neg, individuals_null, col_pos, col_neg, col_null)
savefig(scatter_plot, "scatter_plot.pdf")
savefig(hist_plot, "hist_plot.pdf")

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
