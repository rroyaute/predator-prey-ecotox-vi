using QuadGK
using CSV
using FastGaussQuadrature
using DataFrames
using KernelDensity
using Random
using Statistics
using Distributions
using Plots

# Probability distribution
pdist(y, Y, sigma) = (1 / sqrt(2π * sigma)) * exp(-(y - Y)^2 / (2sigma))

# Attack rate
attack(y, a, alpha, tau) = alpha * exp(-(y - a)^2 / (2tau^2))

# Integral 1 (equivalent to Int1)
function int_1(R, a, alpha, tau, Y, sigma)
    f(y) = attack(y, a, alpha, tau) * pdist(y, Y, sigma)
    val, _ = quadgk(f, -Inf, Inf) # This evaluates the integral of the function f(y) from -Inf to +Inf
    return val
end

# --- LOADING DOSE RESPONSE DATA ---
datadir = joinpath(@__DIR__, "..", "outputs", "data")
mkpath(datadir)
endings = ["pos", "neg", "null"]
filenames = [joinpath(datadir,"dose-response_" * ending * ".csv") for ending in endings]
filename_pos, filename_neg, filename_null = filenames
summary_pos = CSV.read(filename_pos, DataFrame)
summary_neg = CSV.read(filename_neg, DataFrame)
summary_null = CSV.read(filename_null, DataFrame)

function lognormal_params(ybar, sigy)
    sigma2 = log(1 + (sigy / ybar)^2)
    mu = log(ybar) - sigma2 / 2
    return mu, sqrt(sigma2)
end

function propagate_lognormal_gh(f, ybar, sigy; nodes, weights, params_gen)
    mu, sigma = lognormal_params(ybar, sigy)

    ys = exp.(mu .+ sqrt(2) .* sigma .* nodes)
    fy = f.(ys, params_gen.a, params_gen.alpha, params_gen.tau)

    fmean = sum(weights .* fy) / sqrt(pi)
    fvar  = sum(weights .* (fy .- fmean).^2) / sqrt(pi)
    return fmean, sqrt(fvar)
end

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
    a       = 10, # Optimum of alpha is for y=10 now
    h       = 2.5,
    sigma   = sqrt(3)
)

# Warming up functions
_ = attack(params_gen.Y, params_gen.a, params_gen.alpha, params_gen.tau)
nodes, weights = gausshermite(20)

results_pos = [propagate_lognormal_gh(attack, row.yhat_mean, row.yhat_SD; nodes, weights, params_gen) for row in eachrow(summary_pos)]
summary_pos.alpha_mean = first.(results_pos)
summary_pos.alpha_SD = last.(results_pos)

results_neg = [propagate_lognormal_gh(attack, row.yhat_mean, row.yhat_SD; nodes, weights, params_gen) for row in eachrow(summary_neg)]
summary_neg.alpha_mean = first.(results_neg)
summary_neg.alpha_SD = last.(results_neg)

results_null = [propagate_lognormal_gh(attack, row.yhat_mean, row.yhat_SD; nodes, weights, params_gen) for row in eachrow(summary_null)]
summary_null.alpha_mean = first.(results_null)
summary_null.alpha_SD = last.(results_null)

function plot_fy_density(df, xvals, f; n=100_000, params_gen)
    plot()
    for xval in xvals
        row = df[findfirst(==(xval), df.Dose), :]
        mu, sigma = lognormal_params(row.yhat_mean, row.yhat_SD)
        ys = rand(LogNormal(mu, sigma), n)
        fys = f.(ys, params_gen.a, params_gen.alpha, params_gen.tau)

        kd = kde(fys)
        plot!(kd.x, kd.density, label="x = $xval", lw=2)
    end
    xlabel!("y")
    ylabel!("density")
end

plt_density = plot_fy_density(summary_pos, [0.0, 15.0, 30.0], attack; params_gen)
plotdir = joinpath(@__DIR__, "..", "outputs", "figs")
savefig(plt_density, joinpath(plotdir, "test_density_alpha.pdf"))