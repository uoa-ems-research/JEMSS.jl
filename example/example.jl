using JEMSS

println("Loading sim.")
sim = initSim(joinpath(@__DIR__, "input/sim_config.xml"))

println("Simulating.")
simulate!(sim)

println("Writing stats files to \"$(abspath(sim.outputPath))\"")
writeStatsFiles(sim)

# read some stats from file
println("Some statistics from stats_dict.csv output file:")
statsDict = readStatsDictFile(sim.outputFiles["statsDict"].path)
@show statsDict["ambs_avgDailyBusyDurationHours"]
@show statsDict["calls_all_fracResponsesInTime"]

# tables = readTablesFromFile(sim.outputFiles["ambulancesStats"].path)
# sim.outputFiles["callsStats"].path

# # Animate, if you like:
# animate!(sim)

nothing
