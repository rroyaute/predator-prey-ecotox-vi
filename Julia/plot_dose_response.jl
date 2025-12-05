using DataFrames
using Plots
using StatsPlots
using Colors
using LaTeXStrings
using CSV

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

color_pos = RGB(0.01, 0.39, 0.57)
color_neg = RGB(0.78, 0.18, 0.16)
color_null = RGB(0.13, 0.59, 0.12)


plotdir = "../outputs/figs"

# --- LOADING DOSE RESPONSE DATA ---
datadir = "../outputs/data"
endings = ["pos", "neg", "null"]
filenames = [joinpath(datadir,"dose-response_" * ending * ".csv") for ending in endings]
filename_pos, filename_neg, filename_null = filenames
summary_pos = CSV.read(filename_pos, DataFrame)
summary_neg = CSV.read(filename_neg, DataFrame)
summary_null = CSV.read(filename_null, DataFrame)

# --- LOADING STABILITY DATA ---
plotname = "fig_3_repro_julia.pdf"
filenames = ["fig_3_data_julia.csv", "fig_3_d_vals.csv", "fig_3_sigma_vals.csv"]
dataname, dname, sigmaname = filenames
datafile = joinpath(datadir, dataname)
dfile = joinpath(datadir, dname)
sigmafile = joinpath(datadir, sigmaname)
plotfile = joinpath(plotdir, plotname)
df_data = CSV.read(datafile, DataFrame)
d_vals = CSV.read(dfile, DataFrame)[:, 1]
sigma_vals = CSV.read(sigmafile, DataFrame)[:, 1]
data = Matrix(df_data)
println("CSV files loaded successfully.")


# --- GETTING CROPPED VALUES ---
coords_pos = summary_pos[summary_pos.d .<= 2.5, ["yhat_SD", "d"]]
coords_neg = summary_neg[summary_neg.d .<= 2.5, ["yhat_SD", "d"]]
coords_null = summary_null[summary_null.d .<= 2.5, ["yhat_SD", "d"]]

# --- GETTING SPECIAL VALUES ---
x_0 = 4.5
NEC = 1.25
function get_special_values(summary, NEC, x_0)
    coords_NEC = summary[summary.Dose .== NEC, ["yhat_SD", "d"]]
    coords_EC50 = first(summary[summary.d .>= x_0/2, ["yhat_SD", "d"]], 1)
    return coords_NEC, coords_EC50
end
coords_NEC_pos, coords_EC50_pos = get_special_values(summary_pos, NEC, x_0)
coords_NEC_neg, coords_EC50_neg = get_special_values(summary_neg, NEC, x_0)
coords_NEC_null, coords_EC50_null = get_special_values(summary_null, NEC, x_0)
coords_NEC = DataFrame(yhat_SD = [coords_NEC_pos.yhat_SD[1], coords_NEC_neg.yhat_SD[1], coords_NEC_null.yhat_SD[1]],
                       d = [coords_NEC_pos.d[1], coords_NEC_neg.d[1], coords_NEC_null.d[1]])
coords_EC50 = DataFrame(yhat_SD = [coords_EC50_pos.yhat_SD[1], coords_EC50_neg.yhat_SD[1], coords_EC50_null.yhat_SD[1]],
                        d = [coords_EC50_pos.d[1], coords_EC50_neg.d[1], coords_EC50_null.d[1]])



# --- COLOR MAPPING ---
palette = [
    colorant"black",       # 0: Noncoexistence
    colorant"white",       # 1: Coexistence: nonoscillatory
    colorant"lightgray",   # 2: Coexistence: damped oscillations
    colorant"gray"         # 3: Coexistence: limit cycles
]

labels = [
    "Noncoexistence",
    "Coexistence: nonoscillatory",
    "Coexistence: damped oscillations",
    "Coexistence: limit cycles"
]

# Convert integer matrix to colors
color_data = [palette[v+1] for v in data]  # Julia indices start at 1

# Plot
plt = heatmap(
              sigma_vals,
              d_vals,
              color_data,
              yflip = false,
              ylims = (0, 2.5),
            #   yticks = (d_tick_indices, d_labels),
            #   xticks = (sigma_tick_indices, sigma_labels),
              ylabel = "Phenotypic mismatch (d²)",
              xlabel = "Individual variation (σ²)",
              aspect_ratio = 1/1.5,
              size = (1000, 600),
              left_margin = 10Plots.mm,
              right_margin = 20Plots.mm,
              c = palette,
              legend = :outerright,
              legendfontsize = 12,
              grid = false,
              framestyle = :box,
)

# --- ADD LINES OF MEAN AND STD WITH VARYING DOSE ---
plot!(plt, coords_pos.yhat_SD, coords_pos.d,
     color = color_pos,
     linestyle = :solid,
     linewidth = 2,
     label = "Positive cov"
)
plot!(plt, coords_neg.yhat_SD, coords_neg.d,
     color = color_neg,
     linestyle = :solid,
     linewidth = 2,
     label = "Negative cov"
)
plot!(plt, coords_null.yhat_SD, coords_null.d,
     color = color_null,
     linestyle = :solid,
     linewidth = 2,
     label = "No cov"
)

scatter!(plt, coords_NEC.yhat_SD, coords_NEC.d,
    markershape = :cross,
    color = :black,
    label = L"\overline{\mathrm{NEC}}"
)
scatter!(plt, coords_EC50.yhat_SD, coords_EC50.d,
    color = :lightgray, markerstrokewidth = 0,
    label = L"\mathrm{EC}_{50}"
)

# --- ADD LEGEND USING DUMMY POINTS ---
# Plot invisible points with the correct color for legend
for i in 1:4
    scatter!(plt, [NaN], [NaN],
             color = palette[i],
             label = labels[i],
             markershape = :rect,
             markersize=20)
end

plotname = joinpath(plotdir, "background_dose_response_julia.pdf")
savefig(plt, plotname)
println("Plot saved as $plotname.")