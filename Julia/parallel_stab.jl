using Distributed
using Dates
using JLD2

d_vals = range(0, stop=2.5, step=0.005)
sigma_vals = range(0.01, stop=6, step=0.01)

# Add worker processes
addprocs(13)  # 14 cores available
@everywhere using QuadGK, Roots, LinearAlgebra
@everywhere include("method_jacob.jl")

# --- WARM-UP ---
println("Warming up functions...")
warmup_d = d_vals[1]
warmup_sigma = sigma_vals[1]
# call classifying once on the master and once on a worker
classifying(warmup_d, warmup_sigma)
@everywhere classifying($warmup_d, $warmup_sigma) # ensures all processes compile

println("Warm-up complete.")

# --- PREPARE PAIRS ---
pairs = [(d, sigma) for d in d_vals, sigma in sigma_vals]
pairs_flat = vec(pairs)

# --- TIMED PARALLEL COMPUTATION ---
println("Starting parallel computation...")
start_time = now()

results_flat = pmap(pair -> classifying(pair[1], pair[2]), pairs_flat)

end_time = now()
elapsed = end_time - start_time
println("Parallel computation complete.")
println("Elapsed time: $(Dates.value(elapsed)/1e3) seconds (~$(Dates.value(elapsed)/1e3/60) minutes)")

# --- RESHAPE RESULTS ---
data = reshape(results_flat, length(d_vals), length(sigma_vals))

# --- SAVE TO FILE ---
@save "classification_data.jld2" data d_vals sigma_vals
println("Data saved to classification_data.jld2")