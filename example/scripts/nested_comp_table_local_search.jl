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

# Run a local search algorithm to optimise a nested compliance table given some objective function (objFn).
# The local search is simple hill climbing, and the neighbourhood of a nested compliance table
# involves swapping or inserting items (station indices) within/into the table.
# The local search stops when a local optimum (as evaluated by simulating) is found.
# Ambulances are assumed to be identical.

using JEMSS
using Random
using Statistics
using StatsFuns
using Base.Iterators: Stateful

global ObjVal = Float

function nestedCompTableLocalSearch!(sim::Simulation, nestedCompTables::Vector{NestedCompTable}; outputFolder::String="")
    isdir(outputFolder) || mkdir(outputFolder)
    @assert(isdir(outputFolder))
    @assert(!isempty(nestedCompTables))

    # parameters
    global solFilename = "$outputFolder/solutions.csv" # final solutions from each local search
    global nestedCompTablesOutputFilename = "$outputFolder/nested_comp_tables.csv" # save final solutions from each search (saved as priority lists)
    global logFilename = "$outputFolder/log.csv" # log search progress to file
    global sense = :max # :min or :max; direction of optimisation for objective function
    global conf = 0.95 # statistical confidence level
    global doPrint = true
    # reps:
    global minNumReps = 2 # min number of reps to simulate each nested compliance table
    global maxNumRepsList = [2, 5] # = vcat(10 * 2 .^ [0:10;], Inf) # can start with a small maximum, and increase each time local optimum found
    global maxNumRepsIncrease = 1 # = 20 # max number of reps to increase by at a time when comparing two nested compliance tables
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

    # keep track of nested compliance tables tried, and their objective values per replication
    global nestedCompTablesObjVals = Dict{NestedCompTable,Vector{ObjVal}}() # nestedCompTablesObjVals[nestedCompTable][repIndex] gives sim replication objective value

    global nestedCompTablesPeriodStatsList = Dict{NestedCompTable,Vector{SimPeriodStats}}() # nestedCompTablesPeriodStatsList[nestedCompTable][repIndex] gives sim replication period stats

    # save the locally optimal nested compliance tables found for each search
    global nestedCompTableSols = Vector{NestedCompTable}() # nestedCompTableSols[i] is solution for nestedCompTables[i]

    global stationCapacities = [station.capacity for station in sim.stations] # stationCapacities[i] gives ambulance holding capacity of station i

    repeatedLocalSearch!(sim, nestedCompTables)
end

function objValLookup(nestedCompTable::NestedCompTable, repIndex::Int)::ObjVal
    global nestedCompTablesObjVals
    if haskey(nestedCompTablesObjVals, nestedCompTable) && length(nestedCompTablesObjVals[nestedCompTable]) >= repIndex
        return nestedCompTablesObjVals[nestedCompTable][repIndex]
    end
    return nullObjVal
end
function objValsLookup(nestedCompTable::NestedCompTable, repIndices::Vector{Int})::Vector{ObjVal}
    return [objValLookup(nestedCompTable, repIndex) for repIndex in repIndices]
end
function objValsLookup(nestedCompTable::NestedCompTable)::Vector{ObjVal}
    global nestedCompTablesObjVals
    return get(nestedCompTablesObjVals, nestedCompTable, ObjVal[])
end

function setNestedCompTableObjVal!(nestedCompTable::NestedCompTable, repIndex::Int, objVal::ObjVal)
    global nestedCompTablesObjVals
    objVals = get!(nestedCompTablesObjVals, nestedCompTable, ObjVal[])
    @assert(length(objVals) == repIndex - 1) # otherwise will overwrite existing value / will have missing values in objVals
    push!(objVals, objVal)
end

function periodStatsListLookup(nestedCompTable::NestedCompTable)::Vector{SimPeriodStats}
    global nestedCompTablesPeriodStatsList
    return get(nestedCompTablesPeriodStatsList, nestedCompTable, SimPeriodStats[])
end

function setNestedCompTablePeriodStats!(nestedCompTable::NestedCompTable, repIndex::Int, period::SimPeriodStats)
    global nestedCompTablesPeriodStatsList
    periods = get!(nestedCompTablesPeriodStatsList, nestedCompTable, SimPeriodStats[])
    @assert(length(periods) == repIndex - 1) # otherwise will overwrite existing value / will have missing values in periods
    push!(periods, period)
end
function setNestedCompTablePeriodStatsList!(nestedCompTable::NestedCompTable, repIndices::UnitRange{Int}, periods::Vector{SimPeriodStats})
    for (repIndex, period) in zip(repIndices, periods)
        setNestedCompTablePeriodStats!(nestedCompTable, repIndex, period)
    end
end

function resetNestedCompTablesPeriodStatsList!()
    global nestedCompTablesPeriodStatsList
    empty!(nestedCompTablesPeriodStatsList)
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
# Mutates: sim.reps, nestedCompTablesObjVals
function simObjVals!(sim::Simulation, nestedCompTable::NestedCompTable, repIndices::UnitRange{Int})::Tuple{Vector{ObjVal},Vector{SimPeriodStats},Int}
    # look-up any saved data
    objValsFound = objValsLookup(nestedCompTable)
    periodsFound = periodStatsListLookup(nestedCompTable)
    n = min(length(objValsFound), length(periodsFound))
    lookupIndices = intersect(1:n, repIndices)
    simRepIndices = setdiff(repIndices, lookupIndices)

    objVals = ObjVal[]
    periods = SimPeriodStats[]
    if !isempty(simRepIndices)
        # need to simulate for simRepIndices
        global parallel, numThreads
        function runRep!(rep::Simulation)
            reset!(rep)
            applyNestedCompTable!(rep, nestedCompTable)
            simulate!(rep)
        end
        reps = sim.reps[simRepIndices]
        if parallel == true
            runParallel!(runRep!, reps...; numThreads=numThreads)
        else
            for rep in reps
                runRep!(rep)
            end
        end
        periods = getRepsPeriodStatsList(reps)
        objVals = objFn.(periods)
        for (i, repIndex) in enumerate(simRepIndices)
            if length(objValsFound) < repIndex
                setNestedCompTableObjVal!(nestedCompTable, repIndex, objVals[i])
            end
            if length(periodsFound) < repIndex
                setNestedCompTablePeriodStats!(nestedCompTable, repIndex, periods[i])
            end
        end
    end

    objVals = vcat(objValsFound[lookupIndices], objVals)
    periods = vcat(periodsFound[lookupIndices], periods)
    @assert(length(objVals) == length(repIndices))
    @assert(length(periods) == length(repIndices))

    return objVals, periods, length(lookupIndices)
end

function applyNestedCompTable!(sim::Simulation, nestedCompTable::NestedCompTable)
    @assert(!sim.used)
    for i = 1:sim.numAmbs
        setAmbStation!(sim, sim.ambulances[i], sim.stations[nestedCompTable[i]])
    end
    initCompTable!(sim, nestedCompTable)
end

function printNestedCompTable(nestedCompTable::NestedCompTable)
    for i = 1:length(nestedCompTable)
        println("item ", i, ": station ", nestedCompTable[i])
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

global solFileHeader = vcat(["search", "numReps", "objVal", "objValHalfWidth"], sort(collect(keys(simStatsEmpty)))) # ... and nestedCompTable
global logFileHeader = vcat(["search", "maxNumRepsListIndex", "maxNumReps", "iter", "move", "i", "j", "usedMove", "numReps", "numLookups", "searchDurationSeconds", "objVal", "objValHalfWidth", "bestObjVal", "objValDiff", "objValDiffHalfWidth"], sort(collect(keys(simStatsEmpty)))) # ... and nestedCompTable
global logFileDict = Dict{String,Any}([s => "" for s in logFileHeader])

# mutates: file
function fileWriteDlmLine!(file::IOStream, fileHeader::Vector{String}, data::Dict{String,T}, nestedCompTable::NestedCompTable) where {T<:Any}
    line = []
    for header in fileHeader
        push!(line, data[header])
    end
    line = vcat(line, nestedCompTable)
    writeDlmLine!(file, line...)
    flush(file)
end

# perform local search, starting at each of the nested compliance tables provided
# mutates: sim
function repeatedLocalSearch!(sim::Simulation, nestedCompTables::Vector{NestedCompTable})
    global doPrint

    # check that sim is using nested compliance table for move-up
    @assert(sim.moveUpData.useMoveUp == true && sim.moveUpData.moveUpModule == compTableModule)

    # open files for writing solution
    global solFilename, logFilename
    solFile = open(solFilename, "w")
    logFile = open(logFilename, "w")

    # write misc data to files
    global warmUpDuration, periodDuration, conf
    numSearches = length(nestedCompTables)
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

        # initial nested compliance table
        nestedCompTable = nestedCompTables[i]
        doPrint && println("Initial nested compliance table $i (of $numSearches):")
        doPrint && printNestedCompTable(nestedCompTable)

        # perform local search until no improvements can be made, then increase max number of reps allowed
        doPrint && println()
        global minNumReps, maxNumRepsList
        numReps = minNumReps # init
        resetNestedCompTablesPeriodStatsList!()
        for j = 1:length(maxNumRepsList)
            maxNumReps = maxNumRepsList[j] < sim.numReps ? Int(maxNumRepsList[j]) : sim.numReps
            logFileDict["maxNumRepsListIndex"] = j
            logFileDict["maxNumReps"] = maxNumReps

            doPrint && println("Search $i (of $numSearches)")
            doPrint && println("Setting maximum number of reps to ", maxNumReps)
            nestedCompTable, numReps = localSearch!(sim, nestedCompTable, numReps, maxNumReps, logFile)
            if maxNumReps == sim.numReps
                break
            end # increasing max number of reps will not help
        end
        doPrint && println("Finished local search $i (of $numSearches); solution:")
        doPrint && printNestedCompTable(nestedCompTable)

        # save best nested compliance table found for this search
        global nestedCompTableSols
        push!(nestedCompTableSols, nestedCompTable)

        # write solution to file
        objVal, objValHalfWidth = calcMeanAndHalfWidth(objValsLookup(nestedCompTable, [1:numReps;]))
        solFileDict = Dict("search" => i, "numReps" => numReps, "objVal" => objVal, "objValHalfWidth" => objValHalfWidth)
        merge!(solFileDict, getPeriodsStatsDict(periodStatsListLookup(nestedCompTable)[1:numReps]))
        fileWriteDlmLine!(solFile, solFileHeader, solFileDict, nestedCompTable)

        # write solutions as priority list file
        global nestedCompTablesOutputFilename
        writePriorityListsFile(nestedCompTablesOutputFilename, nestedCompTableSols)
    end

    # print out all results from each finished local search
    doPrint && println()
    doPrint && println("Nested compliance tables from completed local searches:")
    for i = 1:numSearches
        doPrint && println()
        doPrint && println("Search: ", i)
        nestedCompTable = nestedCompTableSols[i]
        doPrint && printNestedCompTable(nestedCompTable)
        doPrint && println("objective value = ", mean(objValsLookup(nestedCompTable)))
    end

    close(solFile)
    close(logFile)
end

# Local search of nested compliance table, starting at nestedCompTable.
# Move items in nested compliance table, making changes that optimise the
# objective function (objFn) value for the given sense (:min or :max).
# Continue search until no improvement can be made.
# Mutates: sim, logFile
function localSearch!(sim::Simulation, nestedCompTable::NestedCompTable, nestedCompTableNumReps::Int, maxNumReps::Int, logFile::IOStream)::Tuple{NestedCompTable,Int}

    global minNumReps, maxNumRepsIncrease, logFileHeader, logFileDict, sense, doPrint

    # shorthand
    numAmbs = sim.numAmbs
    numStations = sim.numStations

    # track iterations
    iter = 1 # count iterations performed
    startTime = time()
    getSearchDuration() = round(time() - startTime, digits=2)

    # track nested compliance tables tried
    nestedCompTablesTried = Set{NestedCompTable}()

    # keep track of how many reps were simulated for each nested compliance table
    nestedCompTablesRepsDone = Dict{NestedCompTable,Int}()
    getNestedCompTableRepsDone!(nestedCompTable::NestedCompTable) = get!(nestedCompTablesRepsDone, nestedCompTable, 0)
    updateNestedCompTablesRepsDone!(nestedCompTable::NestedCompTable, n::Int) =
        (nestedCompTablesRepsDone[nestedCompTable] = max(n, getNestedCompTableRepsDone!(nestedCompTable)))

    # calculate objective value for starting point
    numReps = nestedCompTableNumReps # do not just use minNumReps, as nestedCompTable may have already been simulated a number of reps
    objVals, periods, numLookups = simObjVals!(sim, nestedCompTable, 1:numReps)
    (objVal, objValHalfWidth) = calcMeanAndHalfWidth(objVals)
    updateNestedCompTablesRepsDone!(nestedCompTable, numReps)
    push!(nestedCompTablesTried, nestedCompTable)
    doPrint && println("starting objective value: ", objVal)

    # starting solution is current best
    bestNestedCompTable = copy(nestedCompTable)
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
    fileWriteDlmLine!(logFile, logFileHeader, logFileDict, nestedCompTable)

    # local search
    doPrint && println("starting local search")
    endPoint = nothing
    iter += 1
    while true
        for i = (numAmbs+numStations):-1:1, j = 1:numAmbs, move = ["insert", "swap"]
            # starting with i > numAmbs will try insert/swap new items into nested compliance table first (instead of only moving existing items within list); should help with speed

            # change nested compliance table, using re-insertion / swapping
            # if i <= numAmbs, move item in nestedCompTable[i] to nestedCompTable[j]
            # if i > numAmbs, try move slot from station i-numAmbs (if num spare slots > 0)
            nestedCompTable = copy(bestNestedCompTable)
            if move == "insert"
                insert!(nestedCompTable, i, j)
            elseif move == "swap"
                swap!(nestedCompTable, i, j)
            end

            newBestFound = false
            if !in(nestedCompTable, nestedCompTablesTried)
                doPrint && print(move, ": i = ", i, ", j = ", j)

                # simulate until there is a significant difference in objective value between nestedCompTable and bestNestedCompTable,
                # or until running out of replications to simulate
                numReps = minNumReps
                meanDiff = halfWidth = 0.0 # init
                objVals = ObjVal[] # for nestedCompTable
                periods = SimPeriodStats[] # for nestedCompTable
                numLookups = 0
                while numReps <= maxNumReps
                    newObjVals, newPeriods, newNumLookups = simObjVals!(sim, nestedCompTable, (length(objVals)+1):numReps)
                    push!(nestedCompTablesTried, nestedCompTable)
                    push!(objVals, newObjVals...)
                    push!(periods, newPeriods...)
                    numLookups += newNumLookups
                    updateNestedCompTablesRepsDone!(nestedCompTable, numReps)
                    if length(bestObjVals) < length(objVals)
                        push!(bestObjVals, simObjVals!(sim, bestNestedCompTable, (length(bestObjVals)+1):numReps)[1]...)
                        updateNestedCompTablesRepsDone!(bestNestedCompTable, numReps)
                    end
                    meanDiff, halfWidth = calcMeanAndHalfWidth(objVals[1:numReps] - bestObjVals[1:numReps])
                    if numReps == maxNumReps
                        break # ran out of reps, will have to guess which nested compliance table is better based on mean
                    elseif (sense == :max && meanDiff < -halfWidth) || (sense == :min && meanDiff > halfWidth)
                        break # bestNestedCompTable is better, stop
                    elseif (sense == :max && meanDiff > halfWidth) || (sense == :min && meanDiff < -halfWidth)
                        # nestedCompTable may be better
                        # need to simulate nestedCompTable up to much as bestNestedCompTable, to avoid possibility of infinite loop
                        newNumReps = getNestedCompTableRepsDone!(bestNestedCompTable)
                        if numReps < newNumReps
                            numReps = min(maxNumReps, newNumReps)
                        elseif numReps == newNumReps
                            break # done, no need to simulate nestedCompTable further
                        else
                            error("numReps > newNumReps, this should not happen.")
                        end
                    else # not sure which nested compliance table is better, simulate more reps
                        newNumReps = numReps * (halfWidth / meanDiff)^2 # since standard error is roughly proportional to 1/sqrt(n) for n samples
                        if isnan(newNumReps)
                            newNumReps = maxNumReps
                        end # meanDiff was zero, just use all reps
                        numReps = ceil(Int, max(numReps + 1, min(newNumReps, numReps + maxNumRepsIncrease, maxNumReps))) # bound newNumReps
                    end
                end

                objVal, objValHalfWidth = calcMeanAndHalfWidth(objVals)
                bestObjVal = mean(bestObjVals)
                # note that length(objVals) <= length(bestObjVals), as the nested compliance tables may have been simulated different amounts (reps)
                doPrint && println("; objective value: ", objVal, " (best: ", bestObjVal, ")")

                usedMove = false
                if (sense == :max && meanDiff > 0) || (sense == :min && meanDiff < 0) # do not compare objVals and bestObjVals (or their means) as lengths may differ
                    # keep change
                    bestNestedCompTable = nestedCompTable
                    bestObjVals = objVals
                    bestObjVal = objVal
                    usedMove = true
                    newBestFound = true

                    resetNestedCompTablesPeriodStatsList!() # remove periods to save on memory; new solution has almost completely new neighbourhood
                    setNestedCompTablePeriodStatsList!(nestedCompTable, 1:length(periods), periods)

                    doPrint && println("made ", move, ": i = ", i, ", j = ", j, "; new objective value: ", bestObjVal)
                end

                merge!(logFileDict, getPeriodsStatsDict(periods))
                merge!(logFileDict, Dict("iter" => iter, "move" => move, "i" => i, "j" => j, "numReps" => numReps, "numLookups" => numLookups, "usedMove" => Int(usedMove),
                    "searchDurationSeconds" => getSearchDuration(), "objVal" => objVal, "objValHalfWidth" => objValHalfWidth, "bestObjVal" => bestObjVal,
                    "objValDiff" => meanDiff, "objValDiffHalfWidth" => halfWidth))
                fileWriteDlmLine!(logFile, logFileHeader, logFileDict, nestedCompTable)

                iter += 1
            end

            if newBestFound || endPoint === nothing
                endPoint = (i, j, move)
            elseif endPoint == (i, j, move)
                doPrint && println("Found local optimum.")
                return bestNestedCompTable, getNestedCompTableRepsDone!(bestNestedCompTable)
            end
        end
    end
end

# for the given nested compliance table, count how many slots are still available for the station
function countStationSpareSlots(nestedCompTable::NestedCompTable, stationIndex::Int)
    global stationCapacities
    stationSpareSlots = stationCapacities[stationIndex]
    for i in nestedCompTable
        if i == stationIndex
            stationSpareSlots -= 1
        end
    end
    return stationSpareSlots
end

# Insert item at index i into index j in nested compliance table.
# If i > numAmbs, insertion is done with station index i-numAmbs.
# Return true if insertion completed, false otherwise.
function insert!(nestedCompTable::NestedCompTable, i::Int, j::Int)
    if i == j
        return false
    end

    # shorthand
    numAmbs = length(nestedCompTable)
    numStations = length(stationCapacities) # global stationCapacities

    if i <= numAmbs
        @assert(i >= 1)
        # re-insert item in nestedCompTable[i] to nestedCompTable[j]
        temp = nestedCompTable[i]
        if i < j
            for k = i:j-1
                nestedCompTable[k] = nestedCompTable[k+1]
            end
        elseif i > j
            for k = i:-1:j+1
                nestedCompTable[k] = nestedCompTable[k-1]
            end
        else
            error()
        end
        nestedCompTable[j] = temp
    else
        i -= numAmbs # so now i = station index
        @assert(1 <= i && i <= numStations)
        stationSpareSlots = countStationSpareSlots(nestedCompTable, i)
        if stationSpareSlots == 0
            return false
        end
        # insert slot from station i to nestedCompTable[j] (if stationSpareSlots > 0)
        temp = nestedCompTable[numAmbs] # will be bumped from end of nestedCompTable
        for k = numAmbs:-1:j+1
            nestedCompTable[k] = nestedCompTable[k-1]
        end
        nestedCompTable[j] = i
    end

    return true
end

# Swap items at index i and index j in nested compliance table.
# If i > numAmbs, swap is done with station index i-numAmbs.
# Return true if swap completed, false otherwise.
function swap!(nestedCompTable::NestedCompTable, i::Int, j::Int)
    if i == j
        return false
    end

    # shorthand
    numAmbs = length(nestedCompTable)
    numStations = length(stationCapacities) # global stationCapacities

    @assert(1 <= i && i <= numAmbs + numStations)
    @assert(1 <= j && j <= numAmbs)

    # perform swap if stations are different
    if i <= numAmbs
        if nestedCompTable[i] == nestedCompTable[j]
            return false
        end
        (nestedCompTable[i], nestedCompTable[j]) = (nestedCompTable[j], nestedCompTable[i])
    else # elseif i <= numAmbs + numStations
        i -= numAmbs # so now i = station index
        @assert(1 <= i && i <= numStations)
        stationSpareSlots = countStationSpareSlots(nestedCompTable, i)
        if i == nestedCompTable[j] || stationSpareSlots == 0
            return false
        end
        # swap i with nestedCompTable[j]
        nestedCompTable[j] = i
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
nestedCompTableRng = Random.MersenneTwister(0) # useful for reproducing results, if using random restarts
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

# generate nested compliance tables
stationCapacities = [station.capacity for station in sim.stations]
nestedCompTables = [makeRandNestedCompTable(sim.numAmbs, sim.numStations; stationCapacities=stationCapacities, rng=nestedCompTableRng) for i = 1:numSearches]

# set sim to use nested compliance table for move up
setMoveUpModule!(sim, compTableModule)

makeRepsRunnable!(sim) # copy changes to sim (stats, nested compliance table) to sim.reps

# run
println("Starting local search(es).")
nestedCompTableLocalSearch!(sim, nestedCompTables, outputFolder=outputFolder)
println("Total runtime: ", round(time() - t, digits=2), " seconds")
