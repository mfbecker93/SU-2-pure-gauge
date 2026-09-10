using ADerrors

const INPUT = joinpath(@__DIR__, "..", "input", "plaquette_beta2.csv")

lines = readlines(INPUT)

plaquette = Float64[
    parse(Float64, split(line, ",")[2])
    for line in lines[2:end]
]

println("Number of measurements: ", length(plaquette))
println("Mean: ", sum(plaquette) / length(plaquette))

obs = ADerrors.uwreal(plaquette, "plaquette")

ADerrors.uwerr(obs)

cf = obs.cfd[1]

println("length(gamm) = ", length(cf.gamm))

for t in 0:min(29, length(cf.gamm)-1)
    rho = cf.gamm[t + 1] / cf.gamm[1]
    println(t, "  ", rho)
end

println("\n--- ADerrors gamma ---")

for t in 0:9
    println(
        t,
        "  gamma = ",
        obs.cfd[1].gamm[t+1]
    )
end




println("ADerrors mean   = ", ADerrors.value(obs))
println("ADerrors error  = ", ADerrors.err(obs))
println("ADerrors taui   = ", ADerrors.taui(obs, "plaquette"))
println("ADerrors dtaui  = ", ADerrors.dtaui(obs, "plaquette"))
println("ADerrors window = ", ADerrors.window(obs, "plaquette"))

using Printf

mean_val   = ADerrors.value(obs)
error_val  = ADerrors.err(obs)
taui_val   = ADerrors.taui(obs, "plaquette")
dtaui_val  = ADerrors.dtaui(obs, "plaquette")
window_val = ADerrors.window(obs, "plaquette")

println("\\begin{table}[ht]")
println("\\centering")
println("\\begin{tabular}{lc}")
println("\\hline")
println("Parameter & Wert \\\\")
println("\\hline")

@printf("ADerrors mean & %.15f \\\\\n", mean_val)
@printf("ADerrors error & %.15e \\\\\n", error_val)
@printf("ADerrors \$\\tau_i\$ & %.15f \\\\\n", taui_val)
@printf("ADerrors \$\\Delta\\tau_i\$ & %.15f \\\\\n", dtaui_val)
@printf("ADerrors window & %d \\\\\n", window_val)

println("\\hline")
println("\\end{tabular}")
println("\\caption{ADerrors results}")
println("\\label{tab:aderrors}")
println("\\end{table}")




println()
println("--- ADerrors workspace ---")

println("number of fluctuations = ", length(ADerrors.wsg.fluc))
println("map_nob = ", ADerrors.wsg.map_nob)
println("id2str = ", ADerrors.wsg.id2str)
println("str2id = ", ADerrors.wsg.str2id)


