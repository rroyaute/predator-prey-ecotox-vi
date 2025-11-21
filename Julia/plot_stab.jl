
# --- PLOTTING ---
using Plots
using Colors
using JLD2

# --- LOAD DATA ---
@load "classification_data.jld2" data d_vals sigma_vals
println("Data loaded successfully.")

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


# --- TICKS ---
d_step = d_vals[2] - d_vals[1]
sigma_step = sigma_vals[2] - sigma_vals[1]
d_positions = collect(0:0.005:2.5)
sigma_positions = collect(0:0.01:6)
d_tick_step = round(d_step * length(d_positions) / 2.5)
sigma_tick_step = round(sigma_step * length(sigma_positions) / 6)
d_tick_indices = 1:100:length(d_positions)
sigma_tick_indices = 1:100:length(sigma_positions)
d_labels = round.(d_positions[d_tick_indices]; digits=2)
sigma_labels = round.(sigma_positions[sigma_tick_indices]; digits=2)

# Plot
plt = heatmap(
              color_data,
              yflip = false,
              yticks = (d_tick_indices, d_labels),
              xticks = (sigma_tick_indices, sigma_labels),
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

# --- ADD LEGEND USING DUMMY POINTS ---
# Plot invisible points with the correct color for legend
for i in 1:4
    scatter!(plt, [NaN], [NaN],
             color = palette[i],
             label = labels[i],
             markershape = :rect,
             markersize=20)
end

# --- SAVE ---
savefig(plt, "classification_plot.pdf")
println("Plot saved to classification_plot.pdf")