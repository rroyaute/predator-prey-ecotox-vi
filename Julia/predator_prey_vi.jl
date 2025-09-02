using DifferentialEquations
using Plots

function MRmodel!(du, u, p, t)
    x, y = u
    r, K, a, h, e, m = p
    du[1] = dx = r * (1 - x/K) - a * x * y / (1 + a * h * x)
    du[2] = dy = e * a * x * y / (1 + a * h * x) - m * y
end

u0 = [1.0; 0.05] 
p = [0.3; 1.0; 1.76; 1.11; 0.5; 0.1]
tspan = (0.0, 1000.0)
prob = ODEProblem(MRmodel!, u0, tspan, p)
sol = solve(prob)

plot(sol, idxs = (1, 2))
plot(sol, idxs = (0, 2))

prob = DE.ODEProblem(lorenz!, u0, tspan)
sol = DE.solve(prob)