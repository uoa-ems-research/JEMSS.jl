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

# Run a local search algorithm to optimise a priority list given some objective function (objFn).
# The local search is simple hill climbing, and the neighbourhood of a priority list
# involves swapping or inserting items (station indices) within/into the list.
# The local search stops when a local optimum (as evaluated by simulating) is found.
# Ambulances are assumed to be identical.

using JEMSS
using Random
using Statistics
using StatsFuns
using Base.Iterators: Stateful

global ObjVal = Float

function priorityListLocalSearch!(sim::Simulation, priorityLists::Vector{PriorityList}; outputFolder::String="")
    isdir(outputFolder) || mkdir(outputFolder)
    @assert(isdir(outputFolder))
    @assert(!isempty(priorityLists))
    @assert(all(priorityList -> checkPriorityList(priorityList, sim), priorityLists))

    # parameters
    global solFilename = "$outputFolder/solutions.csv" # final solutions from each local search
    global priorityListsOutputFilename = "$outputFolder/priority_lists.csv" # save final priority lists from each search
    global logFilename = "$outputFolder/log.csv" # log search progress to file
    global sense = :max # :min or :max; direction of optimisation for objective function
    global conf = 0.95 # statistical confidence level
    global doPrint = true
    # reps:
    global minNumReps = 2 # min number of reps to simulate each priority list
    global maxNumRepsList = [2, 5] # = vcat(10 * 2 .^ [0:10;], Inf) # can start with a small maximum, and increase each time local optimum found
    global maxNumRepsIncrease = 1 # = 20 # max number of reps to increase by at a time when comparing two priority lists
    # parallel (multi-threading):
    global parallel = false
    global numThreads = 1

    # some parameter checks
    @assert(sense == :min || sense == :max)
    @assert(0 <= conf < 1) # should be between 0 and 1, but not equal to 1 otherwise confidence interval is infinite
    @assert(minNumReps >= 2) # cannot calculate confidence interval for less than 2 reps
    @assert(minNumReps <= sim.numReps)
    @assert(all(maxNumRepsList .>= minNumReps))
    @assert(issorted(maxNumRepsList, lt=<=)) # maxNumRepsList should be strictly increasing
    @assert(maxNumRepsIncrease >= 1)
    @assert(isa(parallel, Bool))
    @assert(numThreads >= 1)

    global nullObjVal = sense == :max ? -Inf : Inf

    # keep track of priority lists tried, and their objective values per replication
    global priorityListsObjVals = Dict{PriorityList,Vector{ObjVal}}() # priorityListsObjVals[priorityList][repIndex] gives sim replication objective value

    global priorityListsPeriodStatsList = Dict{PriorityList,Vector{SimPeriodStats}}() # priorityListsPeriodStatsList[priorityList][repIndex] gives sim replication period stats

    # save the locally optimal priority lists found for each search
    global priorityListSols = Vector{PriorityList}() # priorityListSols[i] is solution for priorityLists[i]

    global stationCapacities = [station.capacity for station in sim.stations] # stationCapacities[i] gives ambulance holding capacity of station i

    repeatedLocalSearch!(sim, priorityLists)
end

function objValLookup(priorityList::PriorityList, repIndex::Int)::ObjVal
    global priorityListsObjVals
    if haskey(priorityListsObjVals, priorityList) && length(priorityListsObjVals[priorityList]) >= repIndex
        return priorityListsObjVals[priorityList][repIndex]
    end
    return nullObjVal
end
function objValsLookup(priorityList::PriorityList, repIndices::Vector{Int})::Vector{ObjVal}
    return [objValLookup(priorityList, repIndex) for repIndex in repIndices]
end
function objValsLookup(priorityList::PriorityList)::Vector{ObjVal}
    global priorityListsObjVals
    return get(priorityListsObjVals, priorityList, ObjVal[])
end

function setPriorityListObjVal!(priorityList::PriorityList, repIndex::Int, objVal::ObjVal)
    global priorityListsObjVals
    objVals = get!(priorityListsObjVals, priorityList, ObjVal[])
    @assert(length(objVals) == repIndex - 1) # otherwise will overwrite existing value / will have missing values in objVals
    push!(objVals, objVal)
end

function periodStatsListLookup(priorityList::PriorityList)::Vector{SimPeriodStats}
    global priorityListsPeriodStatsList
    return get(priorityListsPeriodStatsList, priorityList, SimPeriodStats[])
end

function setPriorityListPeriodStats!(priorityList::PriorityList, repIndex::Int, period::SimPeriodStats)
    global priorityListsPeriodStatsList
    periods = get!(priorityListsPeriodStatsList, priorityList, SimPeriodStats[])
    @assert(length(periods) == repIndex - 1) # otherwise will overwrite existing value / will have missing values in periods
    push!(periods, period)
end
function setPriorityListPeriodStatsList!(priorityList::PriorityList, repIndices::UnitRange{Int}, periods::Vector{SimPeriodStats})
    for (repIndex, period) in zip(repIndices, periods)
        setPriorityListPeriodStats!(priorityList, repIndex, period)
    end
end

function resetPriorityListsPeriodStatsList!()
    global priorityListsPeriodStatsList
    empty!(priorityListsPeriodStatsList)
end

# Return mean and half-width of mean of x.
# Assumes that values in x are from population with normal distribution with unknown standard deviation.
function calcMeanAndHalfWidth(x::Vector{T}; conf::Float=conf)::Tuple{Float,Float} where {T<:Real}
    x = convert(Vector{Float}, x)
    return mean(x), tDistrHalfWidth(x; conf=conf)
end

# Objective value for a single period
function objFn(period::SimPeriodStats)::Float
    return period.call.numResponsesInTime / period.call.numCalls # fraction of calls reached in time
end

# Return objective values from multiple replications (repIndices).
# Uses look-up if possible.
# Mutates: sim.reps, priorityListsObjVals
function simObjVals!(sim::Simulation, priorityList::PriorityList, repIndices::UnitRange{Int})::Tuple{Vector{ObjVal},Vector{SimPeriodStats},Int}
    # look-up any saved data
    objValsFound = objValsLookup(priorityList)
    periodsFound = periodStatsListLookup(priorityList)
    n = min(length(objValsFound), length(periodsFound))
    lookupIndices = intersect(1:n, repIndices)
    simRepIndices = setdiff(repIndices, lookupIndices)

    objVals = ObjVal[]
    periods = SimPeriodStats[]
    if !isempty(simRepIndices)
        # need to simulate for simRepIndices
        global parallel, numThreads
        reps = sim.reps[simRepIndices]
        resetReps!(reps; parallel=parallel, numThreads=numThreads)
        for rep in reps
            applyPriorityList!(rep, priorityList)
        end
        simulateReps!(reps; parallel=parallel, numThreads=numThreads)
        periods = getRepsPeriodStatsList(reps)
        objVals = [objFn(period) for period in periods]
        for (i, repIndex) in enumerate(simRepIndices)
            if length(objValsFound) < repIndex
                setPriorityListObjVal!(priorityList, repIndex, objVals[i])
            end
            if length(periodsFound) < repIndex
                setPriorityListPeriodStats!(priorityList, repIndex, periods[i])
            end
        end
    end

    objVals = vcat(objValsFound[lookupIndices], objVals)
    periods = vcat(periodsFound[lookupIndices], periods)
    @assert(length(objVals) == length(repIndices))
    @assert(length(periods) == length(repIndices))

    return objVals, periods, length(lookupIndices)
end

function applyPriorityList!(sim::Simulation, priorityList::PriorityList)
    @assert(!sim.used)
    for i = 1:sim.numAmbs
        setAmbStation!(sim, sim.ambulances[i], sim.stations[priorityList[i]])
    end
    initPriorityList!(sim, priorityList)
end

function printPriorityList(priorityList::PriorityList)
    for i in eachindex(priorityList)
        println("item ", i, ": station ", priorityList[i])
    end
end

function getPeriodsStatsDict(periods::Vector{SimPeriodStats})::Dict{String,Any}
    global periodDuration, simStatsKeys, simStatsEmpty
    @assert(all(p -> isapprox(p.duration, periodDuration), periods))
    stats = flatten(statsDictFromPeriodStatsList(periods))
    return merge(Dict(["$(key)_mean" => stats[key].mean for key in simStatsKeys]), Dict(["$(key)_halfWidth" => stats[key].halfWidth for key in simStatsKeys]))
end

global simStatsEmpty = flatten(statsDictFromPeriodStatsList(SimPeriodStats[]))
global simStatsKeys = collect(keys(simStatsEmpty))
global simStatsEmpty = merge(Dict(["$(key)_mean" => NaN for key in simStatsKeys]), Dict(["$(key)_halfWidth" => NaN for key in simStatsKeys]))

global solFileHeader = vcat(["search", "numReps", "objVal", "objValHalfWidth"], sort(collect(keys(simStatsEmpty)))) # ... and priorityList
global logFileHeader = vcat(["search", "maxNumRepsListIndex", "maxNumReps", "iter", "move", "i", "j", "usedMove", "numReps", "numLookups", "searchDurationSeconds", "objVal", "objValHalfWidth", "bestObjVal", "objValDiff", "objValDiffHalfWidth"], sort(collect(keys(simStatsEmpty)))) # ... and priorityList
global logFileDict = Dict{String,Any}([s => "" for s in logFileHeader])

# mutates: file
function fileWriteDlmLine!(file::IOStream, fileHeader::Vector{String}, data::Dict{String,T}, priorityList::PriorityList) where {T<:Any}
    line = []
    for header in fileHeader
        push!(line, data[header])
    end
    line = vcat(line, priorityList)
    writeDlmLine!(file, line...)
    flush(file)
end

# perform local search, starting at each of the priority lists provided
# mutates: sim
function repeatedLocalSearch!(sim::Simulation, priorityLists::Vector{PriorityList})
    global doPrint

    # check that sim is using priority list for move-up
    @assert(sim.moveUpData.useMoveUp == true && sim.moveUpData.moveUpModule == priorityListModule)

    # open files for writing solution
    global solFilename, logFilename
    solFile = open(solFilename, "w")
    logFile = open(logFilename, "w")

    # write misc data to files
    global warmUpDuration, periodDuration, conf
    numSearches = length(priorityLists)
    miscTable = Table("miscData", ["numAmbs", "numStations", "numSearches", "conf"];
        rows=[[sim.numAmbs, sim.numStations, numSearches, conf]]) # note that numCalls may vary between sim replications
    repsTable = Table("repsData", ["numReps", "warmUpDuration", "periodDuration"];
        rows=[[sim.numReps, warmUpDuration, periodDuration]])
    writeTablesToFile!(logFile, [miscTable, repsTable])
    writeTablesToFile!(solFile, [miscTable, repsTable])

    # write file headers
    writeDlmLine!(solFile, "solutions")
    writeDlmLine!(logFile, "log")
    writeDlmLine!(solFile, solFileHeader..., ["item_$i stationIndex" for i = 1:sim.numAmbs]...)
    writeDlmLine!(logFile, logFileHeader..., ["item_$i stationIndex" for i = 1:sim.numAmbs]...)
    flush(solFile)
    flush(logFile)

    for i = 1:numSearches
        logFileDict["search"] = i

        # initial priority list
        priorityList = priorityLists[i]
        doPrint && println("Initial priority list $i (of $numSearches):")
        doPrint && printPriorityList(priorityList)

        # perform local search until no improvements can be made, then increase max number of reps allowed
        doPrint && println()
        global minNumReps, maxNumRepsList
        numReps = minNumReps # init
        resetPriorityListsPeriodStatsList!()
        for j in eachindex(maxNumRepsList)
            maxNumReps = maxNumRepsList[j] < sim.numReps ? Int(maxNumRepsList[j]) : sim.numReps
            logFileDict["maxNumRepsListIndex"] = j
            logFileDict["maxNumReps"] = maxNumReps

            doPrint && println("Search $i (of $numSearches)")
            doPrint && println("Setting maximum number of reps to ", maxNumReps)
            priorityList, numReps = localSearch!(sim, priorityList, numReps, maxNumReps, logFile)
            if maxNumReps == sim.numReps
                break
            end # increasing max number of reps will not help
        end
        doPrint && println("Finished local search $i (of $numSearches); solution:")
        doPrint && printPriorityList(priorityList)

        # save best priority list found for this search
        global priorityListSols
        push!(priorityListSols, priorityList)

        # write solution to file
        objVal, objValHalfWidth = calcMeanAndHalfWidth(objValsLookup(priorityList, [1:numReps;]))
        solFileDict = Dict("search" => i, "numReps" => numReps, "objVal" => objVal, "objValHalfWidth" => objValHalfWidth)
        merge!(solFileDict, getPeriodsStatsDict(periodStatsListLookup(priorityList)[1:numReps]))
        fileWriteDlmLine!(solFile, solFileHeader, solFileDict, priorityList)

        global priorityListsOutputFilename
        writePriorityListsFile(priorityListsOutputFilename, priorityListSols)
    end

    # print out all results from each finished local search
    doPrint && println()
    doPrint && println("Priority lists from completed local searches:")
    for i = 1:numSearches
        doPrint && println()
        doPrint && println("Search: ", i)
        priorityList = priorityListSols[i]
        doPrint && printPriorityList(priorityList)
        doPrint && println("objective value = ", mean(objValsLookup(priorityList)))
    end

    close(solFile)
    close(logFile)
end

# Local search of priority list, starting at priorityList.
# Move items in priority list, making changes that optimise the
# objective function (objFn) value for the given sense (:min or :max).
# Continue search until no improvement can be made.
# Mutates: sim, logFile
function localSearch!(sim::Simulation, priorityList::PriorityList, priorityListNumReps::Int, maxNumReps::Int, logFile::IOStream)::Tuple{PriorityList,Int}

    global minNumReps, maxNumRepsIncrease, logFileHeader, logFileDict, sense, doPrint

    # shorthand
    numAmbs = sim.numAmbs
    numStations = sim.numStations

    # track iterations
    iter = 1 # count iterations performed
    startTime = time()
    getSearchDuration() = round(time() - startTime, digits=2)

    # track priority lists tried
    priorityListsTried = Set{PriorityList}()

    # keep track of how many reps were simulated for each priority list
    priorityListsRepsDone = Dict{PriorityList,Int}()
    getPriorityListRepsDone!(priorityList::PriorityList) = get!(priorityListsRepsDone, priorityList, 0)
    updatePriorityListsRepsDone!(priorityList::PriorityList, n::Int) =
        (priorityListsRepsDone[priorityList] = max(n, getPriorityListRepsDone!(priorityList)))

    # calculate objective value for starting point
    numReps = priorityListNumReps # do not just use minNumReps, as priorityList may have already been simulated a number of reps
    objVals, periods, numLookups = simObjVals!(sim, priorityList, 1:numReps)
    (objVal, objValHalfWidth) = calcMeanAndHalfWidth(objVals)
    updatePriorityListsRepsDone!(priorityList, numReps)
    push!(priorityListsTried, priorityList)
    doPrint && println("starting objective value: ", objVal)

    # starting solution is current best
    bestPriorityList = copy(priorityList)
    bestObjVals = objVals
    bestObjVal = objVal

    # write starting point
    for (key, value) in logFileDict
        if !in(key, ("search", "maxNumRepsListIndex", "maxNumReps"))
            logFileDict[key] = "" # reset value
        end
    end
    merge!(logFileDict, getPeriodsStatsDict(periods))
    merge!(logFileDict, Dict("iter" => iter, "numReps" => numReps, "numLookups" => numLookups, "searchDurationSeconds" => getSearchDuration(),
        "objVal" => objVal, "objValHalfWidth" => objValHalfWidth, "bestObjVal" => bestObjVal, "objValDiff" => 0, "objValDiffHalfWidth" => 0))
    fileWriteDlmLine!(logFile, logFileHeader, logFileDict, priorityList)

    # local search
    doPrint && println("starting local search")
    endPoint = nothing
    iter += 1
    while true
        for i = (numAmbs+numStations):-1:1, j = 1:numAmbs, move = ["insert", "swap"]
            # starting with i > numAmbs will try insert/swap new items into priority list first (instead of only moving existing items within list); should help with speed

            # change priority list, using re-insertion / swapping
            # if i <= numAmbs, move item in priorityList[i] to priorityList[j]
            # if i > numAmbs, try move slot from station i-numAmbs (if num spare slots > 0)
            priorityList = copy(bestPriorityList)
            if move == "insert"
                insert!(priorityList, i, j)
            elseif move == "swap"
                swap!(priorityList, i, j)
            end

            newBestFound = false
            if !in(priorityList, priorityListsTried)
                doPrint && print(move, ": i = ", i, ", j = ", j)

                # simulate until there is a significant difference in objective value between priorityList and bestPriorityList,
                # or until running out of replications to simulate
                numReps = minNumReps
                meanDiff = halfWidth = 0.0 # init
                objVals = ObjVal[] # for priorityList
                periods = SimPeriodStats[] # for priorityList
                numLookups = 0
                while numReps <= maxNumReps
                    newObjVals, newPeriods, newNumLookups = simObjVals!(sim, priorityList, (length(objVals)+1):numReps)
                    push!(priorityListsTried, priorityList)
                    push!(objVals, newObjVals...)
                    push!(periods, newPeriods...)
                    numLookups += newNumLookups
                    updatePriorityListsRepsDone!(priorityList, numReps)
                    if length(bestObjVals) < length(objVals)
                        push!(bestObjVals, simObjVals!(sim, bestPriorityList, (length(bestObjVals)+1):numReps)[1]...)
                        updatePriorityListsRepsDone!(bestPriorityList, numReps)
                    end
                    meanDiff, halfWidth = calcMeanAndHalfWidth(objVals[1:numReps] - bestObjVals[1:numReps])
                    if numReps == maxNumReps
                        break # ran out of reps, will have to guess which priority list is better based on mean
                    elseif (sense == :max && meanDiff < -halfWidth) || (sense == :min && meanDiff > halfWidth)
                        break # bestPriorityList is better, stop
                    elseif (sense == :max && meanDiff > halfWidth) || (sense == :min && meanDiff < -halfWidth)
                        # priorityList may be better
                        # need to simulate priorityList up to much as bestPriorityList, to avoid possibility of infinite loop
                        newNumReps = getPriorityListRepsDone!(bestPriorityList)
                        if numReps < newNumReps
                            numReps = min(maxNumReps, newNumReps)
                        elseif numReps == newNumReps
                            break # done, no need to simulate priorityList further
                        else
                            error("numReps > newNumReps, this should not happen.")
                        end
                    else # not sure which priority list is better, simulate more reps
                        newNumReps = numReps * (halfWidth / meanDiff)^2 # since standard error is roughly proportional to 1/sqrt(n) for n samples
                        if isnan(newNumReps)
                            newNumReps = maxNumReps
                        end # meanDiff was zero, just use all reps
                        numReps = ceil(Int, max(numReps + 1, min(newNumReps, numReps + maxNumRepsIncrease, maxNumReps))) # bound newNumReps
                    end
                end

                objVal, objValHalfWidth = calcMeanAndHalfWidth(objVals)
                bestObjVal = mean(bestObjVals)
                # note that length(objVals) <= length(bestObjVals), as the priority lists may have been simulated different amounts (reps)
                doPrint && println("; objective value: ", objVal, " (best: ", bestObjVal, ")")

                usedMove = false
                if (sense == :max && meanDiff > 0) || (sense == :min && meanDiff < 0) # do not compare objVals and bestObjVals (or their means) as lengths may differ
                    # keep change
                    bestPriorityList = priorityList
                    bestObjVals = objVals
                    bestObjVal = objVal
                    usedMove = true
                    newBestFound = true

                    resetPriorityListsPeriodStatsList!() # remove periods to save on memory; new solution has almost completely new neighbourhood
                    setPriorityListPeriodStatsList!(priorityList, 1:length(periods), periods)

                    doPrint && println("made ", move, ": i = ", i, ", j = ", j, "; new objective value: ", bestObjVal)
                end

                merge!(logFileDict, getPeriodsStatsDict(periods))
                merge!(logFileDict, Dict("iter" => iter, "move" => move, "i" => i, "j" => j, "numReps" => numReps, "numLookups" => numLookups, "usedMove" => Int(usedMove),
                    "searchDurationSeconds" => getSearchDuration(), "objVal" => objVal, "objValHalfWidth" => objValHalfWidth, "bestObjVal" => bestObjVal,
                    "objValDiff" => meanDiff, "objValDiffHalfWidth" => halfWidth))
                fileWriteDlmLine!(logFile, logFileHeader, logFileDict, priorityList)

                iter += 1
            end

            if newBestFound || endPoint === nothing
                endPoint = (i, j, move)
            elseif endPoint == (i, j, move)
                doPrint && println("Found local optimum.")
                return bestPriorityList, getPriorityListRepsDone!(bestPriorityList)
            end
        end
    end
end

# for the given priority list, count how many slots are still available for the station
function countStationSpareSlots(priorityList::PriorityList, stationIndex::Int)
    global stationCapacities
    stationSpareSlots = stationCapacities[stationIndex]
    for i in priorityList
        if i == stationIndex
            stationSpareSlots -= 1
        end
    end
    return stationSpareSlots
end

# Insert item at index i into index j in priority list.
# If i > numAmbs, insertion is done with station index i-numAmbs.
# Return true if insertion completed, false otherwise.
function insert!(priorityList::PriorityList, i::Int, j::Int)
    if i == j
        return false
    end

    # shorthand
    numAmbs = length(priorityList)
    numStations = length(stationCapacities) # global stationCapacities

    if i <= numAmbs
        @assert(i >= 1)
        # re-insert item in priorityList[i] to priorityList[j]
        temp = priorityList[i]
        if i < j
            for k = i:j-1
                priorityList[k] = priorityList[k+1]
            end
        elseif i > j
            for k = i:-1:j+1
                priorityList[k] = priorityList[k-1]
            end
        else
            error()
        end
        priorityList[j] = temp
    else
        i -= numAmbs # so now i = station index
        @assert(1 <= i && i <= numStations)
        stationSpareSlots = countStationSpareSlots(priorityList, i)
        if stationSpareSlots == 0
            return false
        end
        # insert slot from station i to priorityList[j] (if stationSpareSlots > 0)
        temp = priorityList[numAmbs] # will be bumped from end of priorityList
        for k = numAmbs:-1:j+1
            priorityList[k] = priorityList[k-1]
        end
        priorityList[j] = i
    end

    return true
end

# Swap items at index i and index j in priority list.
# If i > numAmbs, swap is done with station index i-numAmbs.
# Return true if swap completed, false otherwise.
function swap!(priorityList::PriorityList, i::Int, j::Int)
    if i == j
        return false
    end

    # shorthand
    numAmbs = length(priorityList)
    numStations = length(stationCapacities) # global stationCapacities

    @assert(1 <= i && i <= numAmbs + numStations)
    @assert(1 <= j && j <= numAmbs)

    # perform swap if stations are different
    if i <= numAmbs
        if priorityList[i] == priorityList[j]
            return false
        end
        (priorityList[i], priorityList[j]) = (priorityList[j], priorityList[i])
    else # elseif i <= numAmbs + numStations
        i -= numAmbs # so now i = station index
        @assert(1 <= i && i <= numStations)
        stationSpareSlots = countStationSpareSlots(priorityList, i)
        if i == priorityList[j] || stationSpareSlots == 0
            return false
        end
        # swap i with priorityList[j]
        priorityList[j] = i
    end

    return true
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
priorityListRng = Random.MersenneTwister(0) # useful for reproducing results, if using random restarts
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

# generate priority lists
stationCapacities = [station.capacity for station in sim.stations]
priorityLists = [makeRandPriorityList(sim.numAmbs, sim.numStations; stationCapacities=stationCapacities, rng=priorityListRng) for i = 1:numSearches]

# set sim to use priority list for move up
setMoveUpModule!(sim, priorityListModule)

makeRepsRunnable!(sim) # copy changes to sim (stats, priority list) to sim.reps

# run
println("Starting local search(es).")
priorityListLocalSearch!(sim, priorityLists, outputFolder=outputFolder)
println("Total runtime: ", round(time() - t, digits=2), " seconds")
