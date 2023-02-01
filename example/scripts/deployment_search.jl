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

# Run a search algorithm to optimise ambulance deployment, given some objective function (objFn).
# Given an initial deployment, an ambulance is added where it provides the
# greatest benefit, then an ambulance that provides the least benefit is removed.
# The search stops when no further improvements can be made.
# Ambulances are assumed to be identical, and station capacities are ignored.

using JEMSS
using Random
using Statistics
using StatsFuns
using Base.Iterators: Stateful

const StationsNumAmbs = Vector{Int} # type alias

function deploymentSearch!(sim::Simulation, deployments::Vector{Deployment}; outputFolder::String="")
    isdir(outputFolder) || mkdir(outputFolder)
    @assert(isdir(outputFolder))
    @assert(!isempty(deployments))

    # parameters
    global solFilename = "$outputFolder/solutions.csv" # final solutions from each search
    global deploymentsOutputFilename = "$outputFolder/deployments.csv" # solutions, stored as deployments
    global logFilename = "$outputFolder/log.csv" # log search progress to file
    global sense = :max # :min or :max; direction of optimisation for objective function
    global conf = 0.95 # statistical confidence level
    global doPrint = true
    # parallel (multi-threading):
    global parallel = false
    global numThreads = 1

    # some parameter checks
    @assert(sense == :min || sense == :max)
    @assert(0 <= conf < 1) # should be between 0 and 1, but not equal to 1 otherwise confidence interval is infinite
    @assert(isa(parallel, Bool))
    @assert(numThreads >= 1)

    global nullObjVal = MeanAndHalfWidth(NaN, NaN)

    # keep track of station ambulance counts tried, and their objective values
    global stationsNumAmbsObjVal = Dict{StationsNumAmbs,MeanAndHalfWidth}()
    global stationsNumAmbsStats = Dict{StationsNumAmbs,Dict{String,Any}}()

    # save the station ambulance counts found for each search
    global stationsNumAmbsSols = Vector{StationsNumAmbs}() # stationsNumAmbsSols[i] is solution for deployments[i]

    repeatedSearch!(sim, deployments)
end

function objValLookup(stationsNumAmbs::StationsNumAmbs)::MeanAndHalfWidth
    return get(stationsNumAmbsObjVal, stationsNumAmbs, nullObjVal) # return objective value if found, otherwise nullObjVal
end

function meanAndHalfWidth(x::Vector{T}; conf::Float=conf)::MeanAndHalfWidth where {T<:Real}
    x = convert(Vector{Float}, x)
    return MeanAndHalfWidth(mean(x), tDistrHalfWidth(x; conf=conf))
end

# Objective value for a single period
function objFn(period::SimPeriodStats)::Float
    return period.call.numResponsesInTime / period.call.numCalls # fraction of calls reached in time
end

# For completed simulation replications, calculate and return the objective function value mean and half-width of the mean.
function objFn(sim::Simulation)::MeanAndHalfWidth
    @assert(all(rep -> rep.complete, sim.reps))
    periods = getRepsPeriodStatsList(sim.reps)
    return meanAndHalfWidth([objFn(p) for p in periods])
end

# add ambulance to sim at given station (but not from sim.backup)
function addAmb!(sim::Simulation; stationIndex::Int=1)
    @assert(!sim.used)
    amb = Ambulance() # deepcopy(sim.ambulances[1])
    amb.index = sim.numAmbs + 1
    amb.stationIndex = stationIndex
    JEMSS.initAmbulance!(sim, amb)

    push!(sim.ambulances, amb)
    sim.numAmbs += 1
end

# remove ambulance from sim (but not from sim.backup)
function removeAmb!(sim::Simulation; ambIndex::Int=sim.numAmbs)
    @assert(!sim.used)
    deleteat!(sim.ambulances, ambIndex)
    sim.numAmbs -= 1
    for event in filter(e -> e.ambIndex == ambIndex, sim.eventList)
        JEMSS.deleteEvent!(sim.eventList, event)
    end
end

# apply stationsNumAmbs to sim replications, and return objective value (and half-width) and whether lookup was used
# uses look-up if possible
# mutates: sim, stationsNumAmbsObjVal
function simObjVal!(sim::Simulation, stationsNumAmbs::StationsNumAmbs)::Tuple{Tuple{Float,Float},Bool}
    global stationsNumAmbsObjVal, parallel, numThreads
    objVal = objValLookup(stationsNumAmbs)
    if objVal != nullObjVal
        return ((objVal.mean, objVal.halfWidth), true)
    end
    numAmbsDiff = sum(stationsNumAmbs) - sim.numAmbs # number of ambulances to add to sim
    function runRep!(rep::Simulation)
        resetRep!(rep)
        for i = 1:numAmbsDiff
            addAmb!(rep)
        end
        for i = 1:-numAmbsDiff
            removeAmb!(rep)
        end
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
    objVal = objFn(sim)
    stationsNumAmbsObjVal[stationsNumAmbs] = objVal
    return ((objVal.mean, objVal.halfWidth), false)
end

# print number of ambulances at each station
function printStationsNumAmbs(stationsNumAmbs::StationsNumAmbs)
    for i in eachindex(stationsNumAmbs)
        println("station ", i, ": ", stationsNumAmbs[i], " ambulance(s)")
    end
end

function getSimStats(sim::Simulation)
    # some sim statistics to write to log file
    global periodDuration, simStatsKeys, simStatsEmpty
    if all(r -> r.complete, sim.reps)
        periods = getRepsPeriodStatsList(sim.reps)
        @assert(all(p -> isapprox(p.duration, periodDuration), periods))
        stats = flatten(statsDictFromPeriodStatsList(periods))
        return merge(Dict(["$(key)_mean" => stats[key].mean for key in simStatsKeys]), Dict(["$(key)_halfWidth" => stats[key].halfWidth for key in simStatsKeys]))
    end
    return deepcopy(simStatsEmpty)
end

global simStatsEmpty = flatten(statsDictFromPeriodStatsList(SimPeriodStats[]))
global simStatsKeys = collect(keys(simStatsEmpty))
global simStatsEmpty = merge(Dict(["$(key)_mean" => NaN for key in simStatsKeys]), Dict(["$(key)_halfWidth" => NaN for key in simStatsKeys]))

global solFileHeader = vcat(["search", "objVal", "objValHalfWidth"], sort(collect(keys(simStatsEmpty)))) # ... and stationsNumAmbs
global logFileHeader = vcat(["search", "iter", "oper", "stationIndex", "usedLookup", "searchDurationSeconds", "objVal", "objValHalfWidth", "bestObjVal"], sort(collect(keys(simStatsEmpty)))) # ... and stationsNumAmbs
global logFileDict = Dict{String,Any}([s => "" for s in logFileHeader])

# mutates: file
function fileWriteDlmLine!(file::IOStream, fileHeader::Vector{String}, data::Dict{String,T}, stationsNumAmbs::StationsNumAmbs) where {T<:Any}
    line = []
    for header in fileHeader
        push!(line, data[header])
    end
    line = vcat(line, stationsNumAmbs)
    writeDlmLine!(file, line...)
    flush(file)
end

# perform search, starting at each of the deployments provided
# mutates: sim
function repeatedSearch!(sim::Simulation, deployments::Vector{Deployment})
    global doPrint

    # open files for writing solution
    global solFilename, logFilename
    solFile = open(solFilename, "w")
    logFile = open(logFilename, "w")

    # write misc data to files
    global warmUpDuration, periodDuration, conf
    numSearches = length(deployments)
    miscTable = Table("miscData", ["numAmbs", "numStations", "numSearches", "conf"];
        rows=[[sim.numAmbs, sim.numStations, numSearches, conf]]) # note that numCalls may vary between sim replications
    repsTable = Table("repsData", ["numReps", "warmUpDuration", "periodDuration"];
        rows=[[sim.numReps, warmUpDuration, periodDuration]])
    writeTablesToFile!(logFile, [miscTable, repsTable])
    writeTablesToFile!(solFile, [miscTable, repsTable])

    # write file headers
    writeDlmLine!(solFile, "solutions")
    writeDlmLine!(logFile, "log")
    writeDlmLine!(solFile, solFileHeader..., ["station_$i numAmbs" for i = 1:sim.numStations]...)
    writeDlmLine!(logFile, logFileHeader..., ["station_$i numAmbs" for i = 1:sim.numStations]...)
    flush(solFile)
    flush(logFile)

    for i = 1:numSearches
        logFileDict["search"] = i

        # initial deployment
        stationsNumAmbs = deploymentToStationsNumAmbs(deployments[i], sim.numStations)

        # perform search until no improvements can be made
        doPrint && println()
        doPrint && println("Search ", i, " (of ", numSearches, ")")
        doPrint && printStationsNumAmbs(stationsNumAmbs)
        stationsNumAmbs = search!(sim, stationsNumAmbs, logFile)
        doPrint && println("Finished search $i (of $numSearches); solution:")
        doPrint && printStationsNumAmbs(stationsNumAmbs)

        # save best station ambulance count found for this search
        global stationsNumAmbsSols
        push!(stationsNumAmbsSols, stationsNumAmbs)

        # write solution to file
        objValMeanAndHalfWidth = objValLookup(stationsNumAmbs)
        solFileDict = Dict("search" => i, "objVal" => objValMeanAndHalfWidth.mean, "objValHalfWidth" => objValMeanAndHalfWidth.halfWidth)
        merge!(solFileDict, stationsNumAmbsStats[stationsNumAmbs])
        fileWriteDlmLine!(solFile, solFileHeader, solFileDict, stationsNumAmbs)

        global deploymentsOutputFilename
        writeDeploymentsFile(deploymentsOutputFilename, map(stationsNumAmbsToDeployment, stationsNumAmbsSols), sim.numStations)
    end

    # print out all results from each finished search
    doPrint && println()
    doPrint && println("Station ambulance counts from completed searches:")
    for i = 1:numSearches
        doPrint && println()
        doPrint && println("Search: ", i)
        stationsNumAmbs = stationsNumAmbsSols[i]
        doPrint && printStationsNumAmbs(stationsNumAmbs)
        doPrint && println("objective value = ", objValLookup(stationsNumAmbs).mean)
    end

    close(solFile)
    close(logFile)
end

# Search of ambulance deployment, starting at stationsNumAmbs.
# Given an initial deployment, an ambulance is added where it provides the
# greatest benefit, then an ambulance that provides the least benefit is removed.
# Continue search until no improvement can be made.
# mutates: sim, logFile
function search!(sim::Simulation, stationsNumAmbs::StationsNumAmbs, logFile::IOStream)::StationsNumAmbs

    global logFileHeader, logFileDict, stationsNumAmbsStats, sense, doPrint

    # track iterations
    iter = 1 # count iterations performed
    startTime = time()
    getSearchDuration() = round(time() - startTime, digits=2)

    # calculate objective value for starting point
    (objVal, objValHalfWidth), usedLookup = simObjVal!(sim, stationsNumAmbs)
    doPrint && println("starting objective value: ", objVal)
    bestObjVal = objVal # current best objective value
    bestStationsNumAmbs = copy(stationsNumAmbs)

    # write starting point
    for (key, value) in logFileDict
        if key != "search"
            logFileDict[key] = "" # reset value
        end
    end
    stats = getSimStats(sim)
    if !usedLookup
        stationsNumAmbsStats[stationsNumAmbs] = stats
    end
    merge!(logFileDict, stats)
    merge!(logFileDict, Dict("iter" => iter, "usedLookup" => Int(usedLookup), "searchDurationSeconds" => getSearchDuration(),
        "objVal" => objVal, "objValHalfWidth" => objValHalfWidth, "bestObjVal" => bestObjVal))
    fileWriteDlmLine!(logFile, logFileHeader, logFileDict, stationsNumAmbs)

    # search
    doPrint && println("Starting search.")
    endPoint = nothing
    iter += 1
    numStations = sim.numStations # shorthand
    while true
        bestTempStationsNumAmbs = nothing # will hold new deployment for one more ambulance
        bestTempObjVal = 0
        addStationIndex = 0
        for i = 1:numStations
            # test addition of ambulance to station i
            stationsNumAmbs = copy(bestStationsNumAmbs)
            stationsNumAmbs[i] += 1
            doPrint && print("add: i = $i")

            (objVal, objValHalfWidth), usedLookup = simObjVal!(sim, stationsNumAmbs)
            if !usedLookup
                stats = stationsNumAmbsStats[stationsNumAmbs] = getSimStats(sim)
                merge!(logFileDict, stats)
            else
                merge!(logFileDict, stationsNumAmbsStats[stationsNumAmbs])
                doPrint && print("; done")
            end
            doPrint && println("; objective value: $objVal")

            if i == 1 || (sense == :max && bestTempObjVal < objVal) || (sense == :min && bestTempObjVal > objVal)
                bestTempObjVal = objVal
                bestTempStationsNumAmbs = stationsNumAmbs
                addStationIndex = i
            end

            merge!(logFileDict, Dict("iter" => iter, "oper" => "+", "stationIndex" => i, "usedLookup" => Int(usedLookup),
                "searchDurationSeconds" => getSearchDuration(), "objVal" => objVal, "objValHalfWidth" => objValHalfWidth, "bestObjVal" => bestObjVal))
            fileWriteDlmLine!(logFile, logFileHeader, logFileDict, stationsNumAmbs)
        end
        doPrint && println("added ambulance to station: $addStationIndex")

        rmStationIndex = 0
        for i = 1:numStations
            # test removal of ambulance from station i
            stationsNumAmbs = copy(bestTempStationsNumAmbs)
            stationsNumAmbs[i] -= 1
            if stationsNumAmbs[i] < 0
                continue
            end
            doPrint && print("rm: i = $i")

            (objVal, objValHalfWidth), usedLookup = simObjVal!(sim, stationsNumAmbs)
            if !usedLookup
                stats = stationsNumAmbsStats[stationsNumAmbs] = getSimStats(sim)
                merge!(logFileDict, stats)
            else
                merge!(logFileDict, stationsNumAmbsStats[stationsNumAmbs])
                doPrint && print("; done")
            end
            doPrint && println("; objective value: $objVal")

            if (sense == :max && bestObjVal < objVal) || (sense == :min && bestObjVal > objVal)
                bestObjVal = objVal
                bestStationsNumAmbs = stationsNumAmbs
                rmStationIndex = i
            end

            merge!(logFileDict, Dict("iter" => iter, "oper" => "-", "stationIndex" => i, "usedLookup" => Int(usedLookup),
                "searchDurationSeconds" => getSearchDuration(), "objVal" => objVal, "objValHalfWidth" => objValHalfWidth, "bestObjVal" => bestObjVal))
            fileWriteDlmLine!(logFile, logFileHeader, logFileDict, stationsNumAmbs)
        end

        if rmStationIndex == 0
            # no improvement to deployment, search is complete
            doPrint && println("Finished search.")
            return bestStationsNumAmbs
        else
            doPrint && println("removed ambulance from station: $rmStationIndex")
            iter += 1
        end
    end
end

# mutates: sim
function setSimStatsCapture!(sim::Simulation, periodDurationsIter::Stateful, warmUpDuration::Float)
    stats = sim.stats = deepcopy(sim.backup.stats)
    stats.doCapture = true
    (stats.periodDurationsIter, stats.warmUpDuration) = (periodDurationsIter, warmUpDuration)
    stats.nextCaptureTime = sim.startTime + (stats.warmUpDuration > 0 ? stats.warmUpDuration : first(stats.periodDurationsIter))
    sim.backup.stats = deepcopy(sim.stats)
end

# parameters
configFilename = "sim_config.xml"
outputFolder = "output"
numSearches = 1
deploymentRng = Random.MersenneTwister(0) # useful for reproducing results, if using random restarts
# replication parameters
warmUpDuration = 0.1
periodDuration = 1.0

# parameter checks
@assert(numSearches >= 1)
@assert(warmUpDuration >= 0)
@assert(periodDuration > 0)

# initialise sim
t = time()
@assert(isfile(configFilename))
println("Initialising simulation from config: ", configFilename)
sim = initSim(configFilename, createBackup=true)
@assert(sim.numReps >= 2, "Need at least two simulation replications (sim.reps) for statistics.")

# set sim stats capturing
setSimStatsCapture!(sim, Stateful(periodDuration), warmUpDuration)
periodEndTime = sim.startTime + warmUpDuration + periodDuration
for rep in sim.reps
    # statistics should stop collecting before last call arrives (if not earlier), to avoid cool-down period
    @assert(periodEndTime <= rep.calls[end].arrivalTime)
end

makeRepsRunnable!(sim) # copy changes to sim (stats) to sim.reps

# generate deployments
deployments = makeRandDeployments(sim, numSearches; rng=deploymentRng)

# run
deploymentSearch!(sim, deployments, outputFolder=outputFolder)
println("total runtime: ", round(time() - t, digits=2), " seconds")
