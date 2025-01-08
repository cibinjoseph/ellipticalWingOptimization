using SNOW
using Snopt

include("aero.jl")

function objfunc!(g, x)
    L, D, CL, CD, S = aerodynamics(x)

    g[1] = L - 1.7
    g[2:end] = diff(x)
    return D
end

ns = 25
x0 = ones(ns)
# x0, _ = getEllipticalChords(ns)
lx = 0.0 * ones(ns)
ux = 10.0 * ones(ns)
ng = 1+(ns-1)  # number of constraints
g = zeros(ng)
lg = -Inf * ones(ng)
lg[1] = 0.0
ug = zeros(ng)
options = Options(solver=SNOPT())

xopt, fopt, info = minimize(objfunc!, x0, ng, lx, ux, lg, ug, options)

println("xopt = ", xopt)
println("fopt = ", fopt)
println("info = ", info)

# Verification
if info == :Solve_Succeeded || info == :Solved_To_Acceptable_Level || info[1:8] == "Finished"
    L, D, CL, CD, S = aerodynamics(xopt; output_plot=true)
    println("L = $L")
    println("D = $D")
    println("CL = $CL")
    println("CD = $CD")
end
