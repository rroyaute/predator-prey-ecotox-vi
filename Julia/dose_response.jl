using Distributions
using DataFrames
using LinearAlgebra
using Random
using Statistics
using CSV

"""
Generates the populations and computes the average dose-response from individual dose-response relationships.
"""


# --- DOSE RESPONSE FUNCTIONS ---
function dr_log_rn_fun(dose, ymax, beta, nec) # Threshold dose-response for trait y along dose gradient x
    logyhat = ymax * exp(-exp(beta) * (dose - nec) * (dose > nec)) # y is defined by a maximum trait value y_max and is strictly positive
    return exp(logyhat)
end

function dr_fun_log(dose, ymax, beta, nec)
    if dose < nec
        return log(ymax)
    else
        return log(ymax) - exp(beta) * (dose - nec)
    end
end

dr_fun_log_exp(dose, ymax, beta, nec) = exp(dr_fun_log(dose, ymax, beta, nec))


# --- RANDOM POPULATION GENERATION ---
function generate_individuals(ymax, nec, CVy, CVnec, rho, nID; rng = Random.default_rng())
    mu = [ymax, nec]
    sigmas = [ymax * CVy, nec * CVnec]
    rho_mat = [1.0 rho; rho 1.0]
    sigma_mat = Diagonal(sigmas) * rho_mat * Diagonal(sigmas)
    dist = MvNormal(mu, sigma_mat)
    samples = permutedims(rand(rng, dist, nID))

    return [
        Dict("ymax" => sample[1], "NEC" => sample[2], "ID" => i)
        for (i, sample) in enumerate(eachrow(samples))
    ]
end


# individuals_pos = generate_individuals(y_max, NEC, CV_y, CV_nec, rho, n_id)
# individuals_neg = generate_individuals(y_max, NEC, CV_y, CV_nec, -rho, n_id)
# individuals_null = generate_individuals(y_max, NEC, CV_y, CV_nec, 0, n_id)
# println(individuals_pos)


# --- DOSE RESPONSE SIMULATION FUNCTION ---
function simulate_dose_response(individuals, dose_values, beta, sigma; rng=Random.default_rng())
    records = Vector{Dict{String, Float64}}()

    for ind in individuals
        ymax = ind["ymax"]
        nec = ind["NEC"]
        id = ind["ID"]
        for dose in dose_values
            log_yhat = dr_fun_log(dose, ymax, beta, nec)
            yhat = exp(log_yhat)
            y = rand(rng, LogNormal(log_yhat, sigma))
            push!(records, Dict("ID" => id, "Dose" => float(dose), "log_yhat" => log_yhat, "yhat" => yhat, "y" => y))
        end
    end

    by_dose = Dict{Float64, Vector{Dict{String, Float64}}}()
    for record in records
        push!(get!(by_dose, record["Dose"], Vector{Dict{String, Float64}}()), record)
    end

    summary = Vector{Dict{String, Float64}}()
    for (dose, vals) in by_dose
        yhat_values = [v["yhat"] for v in vals]
        push!(summary, Dict( # This is what is returned: 3 columns (dose, mean and std)
            "Dose" => dose,
            "yhat_mean" => mean(yhat_values),
            "yhat_SD" => std(yhat_values),
        ))
    end

    sort!(summary, by = d -> d["Dose"])
    return DataFrame(summary)
end

# --- PARAMETERS DEFINITION ---
y_max = 10  # trait value for dose = 0 (ymax)
beta = -3.5
Dose = range(0, 100)
NEC = 25
n_id = 10_000
CV_y = 0.1
CV_nec = 0.1
sigma = y_max*0.08
rho = 0.8

Random.seed!(42)


# --- WARMING UP FUNCTIONS ---
println("Warming up functions...")
individuals = generate_individuals(y_max, NEC, CV_y, CV_nec, rho, 10)
summary     = simulate_dose_response(individuals, Dose, beta, sigma)
println("Warm up complete.")

# --- MAIN COMPUTATION ---
println("Starting computation.")
individuals_pos  = generate_individuals(y_max, NEC, CV_y, CV_nec, rho, n_id) # Positive covariance between ymax and NEC
individuals_neg  = generate_individuals(y_max, NEC, CV_y, CV_nec, -rho, n_id) # Negative covariance between ymax and NEC
individuals_null = generate_individuals(y_max, NEC, CV_y, CV_nec, 0, n_id) # ymax and NEC independent

summary_pos  = simulate_dose_response(individuals_pos, Dose, beta, sigma)
summary_neg  = simulate_dose_response(individuals_neg, Dose, beta, sigma)
summary_null = simulate_dose_response(individuals_null, Dose, beta, sigma)

# --- ADDING D VALUES ---
"""
d is either computed as the difference between yhat_mean and y_max which is common to
all 3 populations, or as the difference between yhat_mean and yhat_mean[1] (value of Dose=0),
which is different within the 3 populations (due to random population generation).
"""
function adding_d(summary, y_max)
    summary[!, "d"] = abs.(summary.yhat_mean .- y_max)
    # summary[!, "d"] = abs.(summary.yhat_mean .- summary.yhat_mean[1])
    return summary
end
summary_pos  = adding_d(summary_pos, y_max)
summary_neg  = adding_d(summary_neg, y_max)
summary_null = adding_d(summary_null, y_max)

println("Computation done.")

# --- SAVING TO CSV ---
datadir   = "../outputs/data"
endings   = ["pos", "neg", "null"]
filenames = [joinpath(datadir,"dose-response_" * ending * ".csv") for ending in endings]
filename_pos, filename_neg, filename_null = filenames
CSV.write(filename_pos, summary_pos)
CSV.write(filename_neg, summary_neg)
CSV.write(filename_null, summary_null)
println("CSV files saved:")
println("    -> $filename_pos")
println("    -> $filename_neg")
println("    -> $filename_null")
