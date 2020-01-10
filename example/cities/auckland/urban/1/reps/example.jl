using JEMSS
cd(@__DIR__)

println("Loading sim.")
sim = initSim("input/sim_config.xml")

println("Simulating.")
simulateReps!(sim)

println("Writing stats files.")
writeStatsFiles(sim, sim.reps)
