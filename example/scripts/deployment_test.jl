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

# Read in / generate multiple deployments, and simulate each one.

using JEMSS
using Random

# parameters:
configFilename = "sim_config.xml"
outputFolder = "output"
conf = 0.95 # statistical confidence level, note that sim output files may have different conf value
doPrint = true
# parallel (multi-threading):
parallel = false
numThreads = 1
# additional parameters to use if config does not include deployments:
deploymentsOutputFilename = joinpath(outputFolder, "deployments.csv")
numDeployments = 2 # number of deployments to generate and test
deploymentRng = Random.MersenneTwister(0) # for reproducible results, if generating random deployments

# parameter checks
@assert(isfile(configFilename))
@assert(0 <= conf < 1) # should be between 0 and 1, but not equal to 1 otherwise confidence interval is infinite
@assert(isa(parallel, Bool))
@assert(numThreads >= 1)

isdir(outputFolder) || mkdir(outputFolder)

# load sim
doPrint && println("Loading sim.")
sim = initSim(configFilename, createBackup=true)
(numAmbs, numStations) = (sim.numAmbs, sim.numStations)
@assert(numAmbs > 0 && numStations > 1) # otherwise only one deployment is possible
@assert(sim.numReps >= 2, "Need at least two simulation replications (sim.reps) for statistics.")

# read deployments file
if haskey(sim.inputFiles, "deployments")
    (deployments, numStations2) = readDeploymentsFile(sim.inputFiles["deployments"])
    @assert(numStations == numStations2)
elseif haskey(sim.inputFiles, "deploymentPolicies") # compat
    (deployments, numStations2) = readDeploymentPoliciesFile(sim.inputFiles["deploymentPolicies"])
    @assert(numStations == numStations2)
else
    # generate unique random deployments
    deployments = makeRandDeployments(numAmbs, numStations, numDeployments; rng=deploymentRng)
    doPrint && println("Deployments input file not found in config, generated $numDeployments deployments.")
end
numDeployments = length(deployments)

writeDeploymentsFile(deploymentsOutputFilename, deployments, numStations)

# convert deployments, for ease of use
stationsNumAmbsList = [deploymentToStationsNumAmbs(deployment, numStations) for deployment in deployments]

# check which output files sim will write
simOutputFiles = filter(name -> haskey(sim.outputFiles, name), ["ambulancesStats", "callsStats", "hospitalsStats", "stationsStats", "statsDict"])
simOutputFilepaths = map(name -> sim.outputFiles[name].path, simOutputFiles)

# run simulation for each deployment, save stats
statsPeriodsLists = [] # statsPeriodsLists[i][j] is for deployment i, replication j
statsDicts = []
for (i, stationsNumAmbs) in enumerate(stationsNumAmbsList)
    doPrint && print("\rSimulating deployment $i (of $numDeployments).")
    function runRep!(rep::Simulation)
        resetRep!(rep)
        applyStationsNumAmbs!(rep, stationsNumAmbs)
        simulateRep!(rep)
        @assert(getStationsNumAmbs(rep) == stationsNumAmbs) # check that rep had the current deployment applied
    end
    if parallel == true
        runParallel!(runRep!, sim.reps...; numThreads=numThreads)
    else
        for rep in sim.reps
            runRep!(rep)
        end
    end

    writeStatsFiles(sim, sim.reps)

    # move output files into subfolder
    folder = joinpath(outputFolder, string(i))
    isdir(folder) || mkdir(folder)
    for filepath in simOutputFilepaths
        if isfile(filepath)
            mv(filepath, joinpath(folder, basename(filepath)); force=true)
        end
    end

    # get periods stats
    periods = getRepsPeriodStatsList(sim.reps)
    statsDict = statsDictFromPeriodStatsList(periods, conf=conf)
    push!(statsPeriodsLists, periods)
    push!(statsDicts, statsDict)

    # write deployment to file
    deploymentFilename = "deployment.csv"
    deployment = stationsNumAmbsToDeployment(stationsNumAmbs)
    writeDeploymentsFile(joinpath(folder, deploymentFilename), [deployment], numStations)
end
doPrint && println()
doPrint && println("Done.")

# write objective values (fraction of calls reached in time) to file
objVals = [statsDict["calls"]["all"]["fracResponsesInTime"] for statsDict in statsDicts]
objValsTable = Table("objVals", ["deploymentIndex", "mean", "halfWidth"];
    rows=[[string(i), objVals[i].mean, objVals[i].halfWidth] for i = 1:numDeployments])
writeTablesToFile(joinpath(outputFolder, "obj_vals.csv"), [objValsTable])
