# Step 1: Basic Rosenzweig-MacArthur Predator-Prey Model in Julia
# Let's start with the simplest version and test each component

using DifferentialEquations
using Plots

println("Starting basic Rosenzweig-MacArthur model test...")

# Define the basic Rosenzweig-MacArthur ODE system
# dR/dt = r*R*(1 - R/K) - (a*R*C)/(1 + a*h*R)
# dC/dt = e*(a*R*C)/(1 + a*h*R) - m*C

function rosenzweig_macarthur!(du, u, p, t)
    R, C = u  # Resource and Consumer
    r, K, a, h, e, m = p  # Parameters
    
    # Functional response (Type II)
    functional_response = (a * R) / (1 + a * h * R)
    
    # Resource equation
    du[1] = r * R * (1 - R/K) - functional_response * C
    
    # Consumer equation  
    du[2] = e * functional_response * C - m * C
    
    return nothing
end

# Set basic parameters
params = [
    0.5,   # r: intrinsic growth rate of resource
    1.0,   # K: carrying capacity
    2.0,   # a: attack rate
    0.0,   # h: handling time (start with 0 for Type I response)
    0.5,   # e: conversion efficiency
    0.1    # m: consumer mortality rate
]

# Initial conditions
u0 = [0.8, 0.2]  # [Resource, Consumer]

# Time span
tspan = (0.0, 1000.0)

println("Parameters: r=$(params[1]), K=$(params[2]), a=$(params[3]), h=$(params[4]), e=$(params[5]), m=$(params[6])")
println("Initial conditions: R₀=$(u0[1]), C₀=$(u0[2])")

# Create and solve the ODE problem
prob = ODEProblem(rosenzweig_macarthur!, u0, tspan, params)

println("Solving ODE system...")
sol = solve(prob, Tsit5(), saveat=0.1)

println("Solution status: ", sol.retcode)
println("Number of time points: ", length(sol.t))

# Extract solution for plotting
times = sol.t
resources = [u[1] for u in sol.u]
consumers = [u[2] for u in sol.u]

println("Final values: R = $(round(resources[end], digits=4)), C = $(round(consumers[end], digits=4))")

# Create time series plot
p1 = plot(times, resources, label="Resource (R)", lw=2, color=:green)
plot!(p1, times, consumers, label="Consumer (C)", lw=2, color=:red)
xlabel!(p1, "Time")
ylabel!(p1, "Population")
title!(p1, "Rosenzweig-MacArthur: Time Series")

# Create phase portrait
p2 = plot(resources, consumers, label="Trajectory", lw=2, color=:blue)
scatter!(p2, [u0[1]], [u0[2]], label="Start", color=:green, markersize=6)
scatter!(p2, [resources[end]], [consumers[end]], label="End", color=:red, markersize=6)
xlabel!(p2, "Resource (R)")
ylabel!(p2, "Consumer (C)")
title!(p2, "Phase Portrait")

# Combine plots
combined_plot = plot(p1, p2, layout=(2,1), size=(800, 600))
display(combined_plot)

# Test different parameter sets
println("\n" * "="^50)
println("Testing with Type II functional response (h > 0)...")

# Now test with handling time > 0
params_typeII = [
    0.5,   # r
    1.0,   # K
    2.0,   # a
    0.3,   # h: handling time > 0
    0.5,   # e
    0.1    # m
]

prob2 = ODEProblem(rosenzweig_macarthur!, u0, tspan, params_typeII)
sol2 = solve(prob2, Tsit5(), saveat=0.1)

resources2 = [u[1] for u in sol2.u]
consumers2 = [u[2] for u in sol2.u]

println("Type II solution status: ", sol2.retcode)
println("Type II final values: R = $(round(resources2[end], digits=4)), C = $(round(consumers2[end], digits=4))")

# Compare both solutions
p3 = plot(sol.t, resources, label="Type I: Resource", lw=2, color=:lightgreen, linestyle=:dash)
plot!(p3, sol.t, consumers, label="Type I: Consumer", lw=2, color=:pink, linestyle=:dash)
plot!(p3, sol2.t, resources2, label="Type II: Resource", lw=2, color=:green)
plot!(p3, sol2.t, consumers2, label="Type II: Consumer", lw=2, color=:red)
xlabel!(p3, "Time")
ylabel!(p3, "Population")
title!(p3, "Comparison: Type I vs Type II Functional Response")

display(p3)

println("\nBasic Rosenzweig-MacArthur model test completed successfully!")
println("Ready for next step: adding phenotypic variation...")
