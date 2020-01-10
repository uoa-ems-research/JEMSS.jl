using JEMSS

println("Loading sim.")
sim = initSim(joinpath(@__DIR__, "input/sim_config.xml"))

println("Simulating.")
simulate!(sim)

println("Writing stats files.")
writeStatsFiles(sim)

# # Animate, if you like:
# animate(sim)
