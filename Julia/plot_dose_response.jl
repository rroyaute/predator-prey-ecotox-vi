using DataFrames
using Plots
using StatsPlots
using Colors
using LaTeXStrings

include("dose_response.jl")

function plot_distrib(individuals_pos, individuals_neg, individuals_null,
                     col_pos, col_neg, col_null)

    # Extract {x0, NEC} pairs for scatter plot
    data_pos  = [(ind["x0"], ind["NEC"])  for ind in individuals_pos]
    data_neg  = [(ind["x0"], ind["NEC"])  for ind in individuals_neg]
    data_null = [(ind["x0"], ind["NEC"]) for ind in individuals_null]

    # ---- Scatter plot (ListPlot analogue) ----
    plt_scatter = scatter(
        first.(data_pos),  last.(data_pos),
        color = col_pos, ms = 0.8,
        markerstrokecolor = :transparent,
        markerstrokewidth = 0,
        label = "Positive cov",
        xlabel = "x₀", ylabel = "NEC",
        ratio = :equal,
        framestyle=:box, legend=:topright
    )

    scatter!(
        first.(data_neg), last.(data_neg),
        color = col_neg, ms = 0.8,
        markerstrokecolor = :transparent,
        markerstrokewidth = 0,
        label = "Negative cov"
    )

    scatter!(
        first.(data_null), last.(data_null),
        color = col_null, ms = 0.8,
        markerstrokecolor = :transparent,
        markerstrokewidth = 0,
        label = "No cov"
    )

    # ---- SmoothHistogram analogue using density() ----

    # Extract components for KDE
    pos_x0  = [ind["x0"]  for ind in individuals_pos]
    pos_NEC = [ind["NEC"] for ind in individuals_pos]

    neg_x0  = [ind["x0"]  for ind in individuals_neg]
    neg_NEC = [ind["NEC"] for ind in individuals_neg]

    null_x0 = [ind["x0"]  for ind in individuals_null]
    null_NEC = [ind["NEC"] for ind in individuals_null]

    plt_density = plot(
        xlabel="Value", ylabel="Density",
        framestyle=:box, legend=:topright,
        title="Smoothed Distributions"
    )

    # Positive
    density!(plt_density, pos_x0,  color=col_pos,  lw=2, label="Positive cov (x₀)")
    density!(plt_density, pos_NEC, color=col_pos,  lw=2, linestyle=:dash, label="Positive cov (NEC)")

    # Negative
    density!(plt_density, neg_x0,  color=col_neg,  lw=2, label="Negative cov (x₀)")
    density!(plt_density, neg_NEC, color=col_neg,  lw=2, linestyle=:dash, label="Negative cov (NEC)")

    # Null
    density!(plt_density, null_x0,  color=col_null, lw=2, label="Null cov (x₀)")
    density!(plt_density, null_NEC, color=col_null, lw=2, linestyle=:dash, label="Null cov (NEC)")

    return plt_scatter, plt_density
end

function plot_dose_response(summary_pos, summary_neg, summary_null, col_pos, col_neg, col_null)
    Dose = summary_pos.Dose
    p_mean = plot(Dose, summary_pos.yhat_mean,
    linestyle = :dash, linewidth = 2,
    color = col_pos,
    xlabel = "Concentration",
    ylabel = L"\bar{x}",
    label = ""
    )

    plot!(Dose, summary_neg.yhat_mean,
    linestyle = :dash, linewidth = 2,
    color = col_neg,
    label = ""
    )

    plot!(Dose, summary_null.yhat_mean,
    linestyle = :dash, linewidth = 2,
    color = col_null,
    label = ""
    )

    p_std = twinx()

    plot!(p_std, Dose, summary_pos.yhat_SD,
    linestyle = :solid, linewidth = 2,
    color = col_pos,
    ylabel = L"\sigma",
    label = ""
    )

    plot!(p_std, Dose, summary_neg.yhat_SD,
    linestyle = :solid, linewidth = 2,
    color = col_neg,
    label = ""
    )

    plot!(p_std, Dose, summary_null.yhat_SD,
    linestyle = :solid, linewidth = 2,
    color = col_null,
    label = ""
    )

    plot!(Dose*NaN, Dose*NaN, linestyle = :dash, linewidth = 2, color = col_pos, label = L"Positive cov ($\bar{x}$)")
    plot!(Dose*NaN, Dose*NaN, linestyle = :dash, linewidth = 2, color = col_neg, label = L"Negative cov ($\bar{x}$)")
    plot!(Dose*NaN, Dose*NaN, linestyle = :dash, linewidth = 2, color = col_null, label = L"No cov ($\bar{x}$)")
    plot!(Dose*NaN, Dose*NaN, linestyle = :solid, linewidth = 2, color = col_pos, label = L"Positive cov ($\sigma$)")
    plot!(Dose*NaN, Dose*NaN, linestyle = :solid, linewidth = 2, color = col_neg, label = L"Negative cov ($\sigma$)")
    plot!(Dose*NaN, Dose*NaN, linestyle = :solid, linewidth = 2, color = col_null, label = L"No cov ($\sigma$)")

    plot!(legend = :topright, grid = false)

    return p_mean
end

col_pos = RGB(0.01, 0.39, 0.57)
col_neg = RGB(0.78, 0.18, 0.16)
col_null = RGB(0.13, 0.59, 0.12)

x_0 = 4.5
beta = -0.5
Dose = range(0, step = 0.05, length = 101)
NEC = 1.25
n_id = 10_000
CV_x = 1/x_0
CV_nec = 0.2
rho = 0.9
sigma = 0.08

Random.seed!(42)

# --- Warming up of functions --- #

println("Warming up of functions...")

individuals_pos = generate_individuals(x_0, NEC, CV_x, CV_nec, rho, 10)
individuals_neg = generate_individuals(x_0, NEC, CV_x, CV_nec, -rho, 10)
individuals_null = generate_individuals(x_0, NEC, CV_x, CV_nec, 0, 10)

scatter_plot, hist_plot = plot_distrib(individuals_pos, individuals_neg, individuals_null, col_pos, col_neg, col_null)

summary_pos = simulate_dose_response(individuals_pos, Dose, beta, sigma)
summary_neg = simulate_dose_response(individuals_neg, Dose, beta, sigma)
summary_null = simulate_dose_response(individuals_null, Dose, beta, sigma)

summary_plot = plot_dose_response(summary_pos, summary_neg, summary_null, col_pos, col_neg, col_null)

println("Warming up complete.")

# --- Main code --- #

individuals_pos = generate_individuals(x_0, NEC, CV_x, CV_nec, rho, n_id)
individuals_neg = generate_individuals(x_0, NEC, CV_x, CV_nec, -rho, n_id)
individuals_null = generate_individuals(x_0, NEC, CV_x, CV_nec, 0, n_id)

scatter_plot, hist_plot = plot_distrib(individuals_pos, individuals_neg, individuals_null, col_pos, col_neg, col_null)

summary_pos = simulate_dose_response(individuals_pos, Dose, beta, sigma)
summary_neg = simulate_dose_response(individuals_neg, Dose, beta, sigma)
summary_null = simulate_dose_response(individuals_null, Dose, beta, sigma)

summary_plot = plot_dose_response(summary_pos, summary_neg, summary_null, col_pos, col_neg, col_null)


plot_dir = "../outputs/figs"
savefig(scatter_plot, joinpath(plot_dir, "scatter_plot_distrib_x_NEC.pdf"))
println("Scatter plot of distributions saved as scatter_plot_distrib_x_NEC.pdf.")
savefig(hist_plot, joinpath(plot_dir, "hist_plot_distrib_x_NEC.pdf"))
println("Histrogram of x and NEC distributions saved as hist_plot_distrib_x_NEC.pdf.")
savefig(summary_plot, joinpath(plot_dir, "summary_plot_dose_response.pdf"))
println("Dose-response of x (mean and std) saved as summary_plot_dose_response.pdf.")