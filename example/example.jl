using AmbulanceSim

path = @__DIR__

# generate artificial simulation input files
println("\n=== Generating files ===")
genConfigFilename = joinpath(path, "gen_config.xml")
runGenConfig(genConfigFilename)

# create and run simulation using generated files
println("\n=== Simulating with generated files ===")
simConfigFilename = joinpath(path, "sim_config.xml")
sim = initSimulation(simConfigFilename)
simulate!(sim)

# print some basic statistics
println("\n=== Printing some simulation statistics ===")
printSimStats(sim)

# reset simulation and run again (useful for restarting animation)
println("\n=== Re-running simulation ===")
resetSim!(sim)
simulate!(sim)

# create and run simulation again, this time writing output files
println("\n=== Simulating with generated files, writing output ===")
simConfigFilename = joinpath(path, "sim_config.xml")
sim = initSimulation(simConfigFilename; allowWriteOutput = true)
openOutputFiles!(sim)
simulate!(sim)
writeStatsFiles!(sim)
closeOutputFiles!(sim)

# create and run simulation again, resimulating based on previously created events output file
println("\n=== Resimulating based on output/events file ===")
sim = initSimulation(simConfigFilename; allowResim = true)
simulate!(sim)
