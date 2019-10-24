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

@info("Todo: use priorityListsOutputFilename.")

function priorityListLocalSearch!(sim::Simulation; outputFolder::String = "", priorityLists::Vector{PriorityList} = [])

isdir(outputFolder) || mkdir(outputFolder)
@assert(isdir(outputFolder))

@assert(all(priorityList -> checkPriorityList(priorityList, sim), priorityLists))

# parameters
global solFilename = "$outputFolder/solutions.csv" # final solutions from each local search
global priorityListsOutputFilename = "$outputFolder/priority_lists.csv" # save final priority lists from each search
global logFilename = "$outputFolder/log.csv" # log search progress to file
global sense = :max # :min or :max; direction of optimisation for objective function
global conf = 0.95 # statistical confidence level
global doPrint = true

# some parameter checks
@assert(sense == :min || sense == :max)
@assert(0 <= conf < 1) # should be between 0 and 1, but not equal to 1 otherwise confidence interval is infinite

global numSearches
@assert(numSearches == length(priorityLists))

global nullObjVal = MeanAndHalfWidth(NaN,NaN)

# keep track of priority lists tried, and their objective values
global priorityListsObjVal = Dict{PriorityList, MeanAndHalfWidth}()
global priorityListsStats = Dict{PriorityList, Dict{String, Any}}()
function objValLookup(priorityList::PriorityList)::MeanAndHalfWidth
	return get(priorityListsObjVal, priorityList, nullObjVal) # return objective value if found, otherwise nullObjVal
end

function meanAndHalfWidth(x::Vector{T}; conf::Float = conf)::MeanAndHalfWidth where T <: Real
	x = convert(Vector{Float}, x)
	return MeanAndHalfWidth(mean(x), tDistrHalfWidth(x; conf = conf))
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

# run sim replications with priorityList, and return objective value (and half-width) and whether lookup was used
# uses look-up if possible
# mutates: sim, priorityListsObjVal
function simObjVal!(sim::Simulation, priorityList::PriorityList)::Tuple{MeanAndHalfWidth,Bool}
	global priorityListsObjVal
	objVal = objValLookup(priorityList)
	if objVal != nullObjVal return (objVal, true) end
	for rep in sim.reps
		reset!(sim)
		for i = 1:sim.numAmbs
			setAmbStation!(sim, sim.ambulances[i], sim.stations[priorityList[i]])
		end
		initPriorityList!(sim, priorityList)
		simulateRep!(sim, rep)
	end
	objVal = objFn(sim)
	priorityListsObjVal[priorityList] = objVal
	return (objVal, false)
end

# save the locally optimal station ambulance counts found for each search
global priorityListSols = Vector{PriorityList}() # priorityListSols[i] is solution for priorityLists[i]

global stationCapacities = Int[] # stationCapacities[i] gives ambulance holding capacity of station i; will populate after sim is initialised

function printPriorityList(priorityList::PriorityList)
	for i = 1:length(priorityList)
		println("item ", i, ": station ", priorityList[i])
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

global solFileHeader = ["search", "objVal", "objValHalfWidth"] # ... and priorityList
global logFileHeader = vcat(["search", "iter", "move", "i", "j", "usedLookup", "usedMove", "searchDurationSeconds", "objVal", "objValHalfWidth", "bestObjVal"], sort(collect(keys(simStatsEmpty)))) # ... and priorityList
global logFileDict = Dict{String,Any}([s => "" for s in logFileHeader])

# mutates: file
function fileWriteDlmLine!(file::IOStream, fileHeader::Vector{String}, data::Dict{String,T}, priorityList::PriorityList) where T <: Any
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
	global numSearches, warmUpDuration, periodDuration, conf
	miscTable = Table("miscData", ["numAmbs", "numStations", "numSearches", "conf"];
		rows = [[sim.numAmbs, sim.numStations, numSearches, conf]]) # note that numCalls may vary between sim replications
	repsTable = Table("repsData", ["numReps", "warmUpDuration", "periodDuration"];
		rows = [[sim.numReps, warmUpDuration, periodDuration]])
	writeTablesToFile!(logFile, [miscTable, repsTable])
	writeTablesToFile!(solFile, [miscTable, repsTable])
	
	# write file headers
	writeDlmLine!(solFile, "solutions")
	writeDlmLine!(logFile, "log")
	writeDlmLine!(solFile, solFileHeader..., ["item_$i stationIndex" for i = 1:sim.numAmbs]...)
	writeDlmLine!(logFile, logFileHeader..., ["item_$i stationIndex" for i = 1:sim.numAmbs]...)
	flush(solFile)
	flush(logFile)
	
	global stationCapacities
	stationCapacities = [station.capacity for station in sim.stations]
	
	for i = 1:numSearches
		logFileDict["search"] = i
		
		# initial priority list
		priorityList = priorityLists[i]
		
		# perform local search until no improvements can be made
		doPrint && println()
		doPrint && println("Search ", i, " (of ", numSearches, ")")
		doPrint && printPriorityList(priorityList)
		priorityList = localSearch!(sim, priorityList, logFile)
		doPrint && println("Finished local search $i (of $numSearches); solution:")
		doPrint && printPriorityList(priorityList)
		
		# save best priority list found for this search
		global priorityListSols
		push!(priorityListSols, priorityList)
		
		# write solution to file
		objValMeanAndHalfWidth = objValLookup(priorityList)
		solFileDict = Dict("search" => i, "objVal" => objValMeanAndHalfWidth.mean, "objValHalfWidth" => objValMeanAndHalfWidth.halfWidth)
		fileWriteDlmLine!(solFile, solFileHeader, solFileDict, priorityList)
		
		# global priorityListsOutputFilename
		# writePriorityListsFile(priorityListsOutputFilename, priorityListSols)
	end
	
	# print out all results from each finished local search
	doPrint && println()
	doPrint && println("Priority lists from completed local searches:")
	for i = 1:numSearches
		doPrint && println()
		doPrint && println("Search: ", i)
		priorityList = priorityListSols[i]
		doPrint && printPriorityList(priorityList)
		doPrint && println("objective value = ", objValLookup(priorityList).mean)
	end
	
	close(solFile)
	close(logFile)
end

# Local search of priority list, starting at priorityList.
# Move ambulances between stations, making changes that optimise
# the objective function (objFn) for the given sense (:min or :max).
# Continue search until no improvement can be made.
# Mutates: sim, logFile
function localSearch!(sim::Simulation, priorityList::PriorityList, logFile::IOStream)::PriorityList
	
	global logFileHeader, logFileDict, priorityListsStats, sense, doPrint
	
	# track iterations
	iter = 1 # count iterations performed
	startTime = time()
	getSearchDuration() = round(time() - startTime, digits = 2)
	
	# calculate objective value for starting point
	objValMeanAndHalfWidth, usedLookup = simObjVal!(sim, priorityList)
	(objVal, objValHalfWidth) = (objValMeanAndHalfWidth.mean, objValMeanAndHalfWidth.halfWidth)
	doPrint && println("starting objective value: ", objVal)
	bestObjVal = objVal # current best objective value
	bestPriorityList = copy(priorityList)
	
	# write starting point
	for (key, value) in logFileDict
		if key != "search"
			logFileDict[key] = "" # reset value
		end
	end
	stats = getSimStats(sim)
	if !usedLookup priorityListsStats[priorityList] = stats end
	merge!(logFileDict, stats)
	merge!(logFileDict, Dict("iter" => iter, "usedLookup" => Int(usedLookup), "searchDurationSeconds" => getSearchDuration(),
		"objVal" => objVal, "objValHalfWidth" => objValHalfWidth, "bestObjVal" => bestObjVal))
	fileWriteDlmLine!(logFile, logFileHeader, logFileDict, priorityList)
	
	# local search
	doPrint && println("starting local search")
	endPoint = nothing
	iter += 1
	numStations = sim.numStations # shorthand
	numAmbs = sim.numAmbs # shorthand
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
			if priorityList != bestPriorityList
				doPrint && print(move, ": i = ", i, ", j = ", j)
				
				objValMeanAndHalfWidth, usedLookup = simObjVal!(sim, priorityList)
				(objVal, objValHalfWidth) = (objValMeanAndHalfWidth.mean, objValMeanAndHalfWidth.halfWidth)
				if !usedLookup
					stats = priorityListsStats[priorityList] = getSimStats(sim)
					merge!(logFileDict, stats)
				else
					merge!(logFileDict, priorityListsStats[priorityList])
					doPrint && print("; done")
				end
				
				doPrint && println("; objective value: ", objVal, " (best: ", bestObjVal, ")")
				
				usedMove = false
				if (sense == :max && bestObjVal < objVal) || (sense == :min && bestObjVal > objVal)
					bestObjVal = objVal
					# keep change (move ambulance from station i to j)
					bestPriorityList = priorityList
					usedMove = true
					newBestFound = true
					
					doPrint && println("made ", move, ": i = ", i, ", j = ", j, "; new objective value: ", bestObjVal)
				end
				
				merge!(logFileDict, Dict("iter" => iter, "move" => move, "i" => i, "j" => j, "usedLookup" => Int(usedLookup), "usedMove" => Int(usedMove),
					"searchDurationSeconds" => getSearchDuration(), "objVal" => objVal, "objValHalfWidth" => objValHalfWidth, "bestObjVal" => bestObjVal))
				fileWriteDlmLine!(logFile, logFileHeader, logFileDict, priorityList)
				
				iter += 1
			end
			
			if newBestFound || endPoint == nothing
				endPoint = (i,j,move)
			elseif endPoint == (i,j,move)
				doPrint && println("Found local optimum.")
				return bestPriorityList
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
	if i == j return false end
	
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
	if i == j return false end
	
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

repeatedLocalSearch!(sim, priorityLists)

end # priorityListLocalSearch!

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
sim = initSim(configFilename, createBackup = true)
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
priorityLists = [makeRandPriorityList(sim.numAmbs, sim.numStations; stationCapacities = stationCapacities, rng = priorityListRng) for i = 1:numSearches]

# set sim to use priority list for move up
for tempSim in (sim, sim.backup)
	moveUpData = getfield(tempSim, :moveUpData)
	for (fname, val) in ((:useMoveUp, true), (:moveUpModule, priorityListModule))
		setfield!(moveUpData, fname, val)
	end
end

# run
priorityListLocalSearch!(sim, outputFolder = outputFolder, priorityLists = priorityLists)
println("total runtime: ", round(time()-t, digits = 2), " seconds")
