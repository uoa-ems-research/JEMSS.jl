using JEMSS

println("Loading sim.")
sim = initSim(joinpath(@__DIR__, "input/sim_config.xml"))

println("Simulating.")
simulateReps!(sim)

println("Writing stats files.")
writeStatsFiles(sim, sim.reps)
