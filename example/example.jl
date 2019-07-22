##########################################################################
# Copyright 2017 Samuel Ridler.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##########################################################################

using JEMSS

path = @__DIR__

# generate artificial simulation input files
println("\n=== Generating files ===")
genConfigFilename = joinpath(path, "gen_config.xml")
runGenConfig(genConfigFilename; overwriteOutputPath = true)

# create and run simulation using generated files
println("\n=== Simulating with generated files ===")
simConfigFilename = joinpath(path, "sim_config.xml")
sim = initSim(simConfigFilename, doPrint = false)
simulate!(sim)

# print some basic statistics
# println("\n=== Printing some simulation statistics ===")
println()
printSimStats(sim)

# reset simulation and run again (useful for restarting animation)
println("\n=== Re-running simulation ===")
reset!(sim)
simulate!(sim)

# create and run simulation again, this time writing output files
println("\n=== Simulating with generated files, writing output ===")
simConfigFilename = joinpath(path, "sim_config.xml")
sim = initSim(simConfigFilename; allowWriteOutput = true, doPrint = false)
openOutputFiles!(sim)
simulate!(sim)
writeOutputFiles(sim)
closeOutputFiles!(sim)

# create and run simulation again, resimulating based on previously created events output file
println("\n=== Resimulating based on output/events file ===")
sim = initSim(simConfigFilename; allowResim = true, doPrint = false)
simulate!(sim)

nothing # return value
