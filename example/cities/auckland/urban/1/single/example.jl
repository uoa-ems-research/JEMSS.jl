using JEMSS
cd(@__DIR__)

println("Loading sim.")
sim = initSim("input/sim_config.xml")

println("Simulating.")
simulate!(sim)

println("Writing stats files.")
writeStatsFiles(sim)

# # Animate, if you like:
# animate(sim)
