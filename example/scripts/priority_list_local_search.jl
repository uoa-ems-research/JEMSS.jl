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
# Ambulances are assumed to be identical, and station capacities are ignored.

using JEMSS
using Random
using Statistics

# parameters:
const configFilename = "sim_config.xml"
const outputFolder = "output"
const solFilename = "$outputFolder/solutions.csv" # final solutions from each local search
const logFilename = "$outputFolder/log.csv" # log search progress to file
const ObjVal = Int # type alias, return type of objFn
const nullObjVal = -1 # depends on objective function, see objFn
const sense = :max # :min or :max; direction of optimisation for objective function
priorityLists = [] # leave empty (and set numSearches) if generating random priority lists for random restarts
const numSearches = isempty(priorityLists) ? 1 : length(priorityLists) # number of local searches to perform
const priorityListRng = Random.MersenneTwister(0) # useful for reproducing results, if using random restarts

# some parameter checks
@assert(isfile(configFilename))
@assert(isdir(outputFolder))
@assert(ObjVal <: Real)
@assert(typeof(nullObjVal) <: ObjVal)
@assert(sense == :min || sense == :max)
@assert((isempty(priorityLists) && numSearches >= 1) || (length(priorityLists) == numSearches))
@assert(all(priorityList -> checkPriorityList(priorityList), priorityLists))

stationCapacities = [] # stationCapacities[i] gives ambulance holding capacity of station i; will populate after sim is initialised
priorityListSols = Vector{PriorityList}() # solution of each complete local search

# keep track of priority lists tried, and their objective values
priorityListsObjVal = Dict{PriorityList, ObjVal}()
function objValLookup(priorityList::PriorityList)::ObjVal
	return get(priorityListsObjVal, priorityList, nullObjVal) # return objective value if found, otherwise nullObjVal
end

# for a completed simulation, calculate and return the objective function value
function objFn(sim::Simulation)::ObjVal
	@assert(sim.complete)
	return countCallsReachedInTime(sim) # countCallsReachedInTime is from JEMSS
end

# get objective value for sim, applying priorityList
# only use this if priorityListsObjVal does not already contain priorityList
# mutates: sim, priorityListsObjVal
function simObjVal!(sim::Simulation, priorityList::PriorityList)::ObjVal
	@assert(!haskey(priorityListsObjVal, priorityList))
	reset!(sim)
	for i = 1:sim.numAmbs
		setAmbStation!(sim.ambulances[i], sim.stations[priorityList[i]])
	end
	initPriorityList!(sim, priorityList)
	simulateToEnd!(sim)
	objVal = objFn(sim)
	priorityListsObjVal[priorityList] = objVal
	# do not reset sim before returning, as other stats may needed after calling this function
	return objVal
end

function printPriorityList(priorityList::PriorityList)
	for i = 1:length(priorityList)
		println("item ", i, ": station ", priorityList[i])
	end
end

function simStats(sim::Simulation)
	# some sim statistics to write to log file
	stats = Dict{String,Any}()
	stats["totalAmbTravelDuration"] = ""
	stats["totalAmbBusyDuration"] = ""
	stats["avgResponseDurationMinutes"] = ""
	stats["callsReachedInTime"] = ""
	if sim.complete
		stats["totalAmbTravelDuration"] = sum(amb -> amb.totalTravelDuration, sim.ambulances)
		stats["totalAmbBusyDuration"] = sum(amb -> amb.totalBusyDuration, sim.ambulances)
		@assert(all(call -> call.responseDuration != nullTime, sim.calls))
		stats["avgResponseDurationMinutes"] = mean(call -> call.responseDuration, sim.calls) * 24 * 60
		stats["callsReachedInTime"] = sum(call -> call.responseDuration <= sim.targetResponseDurations[Int(call.priority)], sim.calls)
	end
	return stats
end

logFileHeader = ["search", "iter", "move", "i", "j", "usedLookup", "usedMove", "objVal", "bestObjVal", "totalAmbTravelDuration", "totalAmbBusyDuration", "avgResponseDurationMinutes", "callsReachedInTime", "iterTimeSeconds"] # ... and priority list
logFileDict = Dict{String,Any}([s => "" for s in logFileHeader])

function logFileWriteDlmLine!(logFile::IOStream, data::Dict{String,Any}, priorityList::PriorityList)
	line = []
	for header in logFileHeader
		push!(line, data[header])
	end
	line = vcat(line, priorityList)
	writeDlmLine!(logFile, line...)
	flush(logFile)
end

# perform local search, starting at each of the priority lists provided
function repeatedLocalSearch()
	println("Initialising simulation from config: ", configFilename)
	sim = initSim(configFilename; createBackup = false, doPrint = false)
	
	# change move up module index to priority_list
	sim.moveUpData.moveUpModule = priorityListModule
	sim.moveUpData.useMoveUp = true
	
	backup!(sim)
	
	global stationCapacities = [station.capacity for station in sim.stations]
	
	# open files for writing solution
	solFile = open(solFilename, "w")
	logFile = open(logFilename, "w")
	
	# write file headers
	writeDlmLine!(solFile, "search", "objVal", ["item_$i stationIndex" for i = 1:sim.numAmbs]...)
	writeDlmLine!(logFile, logFileHeader..., ["item_$i stationIndex" for i = 1:sim.numAmbs]...)
	flush(solFile)
	flush(logFile)
	
	global priorityLists, priorityListSols
	if isempty(priorityLists)
		println("generating $numSearches random priority list(s)")
		for i = 1:numSearches
			priorityList = makeRandPriorityList(sim.numAmbs, sim.numStations;
				stationCapacities = stationCapacities, rng = priorityListRng)
			push!(priorityLists, priorityList)
		end
	end
	@assert(all(priorityList -> checkPriorityList(priorityList, sim), priorityLists))
	
	for i = 1:numSearches
		logFileDict["search"] = i
		
		# perform local search until no improvements can be made
		println()
		println("Search ", i, " (of ", numSearches, ")")
		priorityList = localSearch!(sim, priorityLists[i], logFile)
		
		# save best ambulance to station allocation found for this search iteration
		push!(priorityListSols, priorityList)
		
		# write solution to file
		writeDlmLine!(solFile, i, objValLookup(priorityList), priorityList...)
		flush(solFile)
	end
	
	# print out all results from each finished local search
	println()
	println("Priority lists from completed local searches:")
	for i = 1:numSearches
		println()
		println("Iteration: ", i)
		priorityList = priorityListSols[i]
		printPriorityList(priorityList)
		println("objective value = ", objValLookup(priorityList))
	end
	
	close(solFile)
	close(logFile)
end

# Local search of priority list.
# Move items in priority list, making changes that optimise
# the objective function (objFn) for the given sense (:min or :max).
# Continue search until no improvement can be made.
# Mutates: sim, logFile
function localSearch!(sim::Simulation, priorityList::PriorityList, logFile::IOStream)
	
	# shorthand
	numAmbs = sim.numAmbs
	numStations = sim.numStations
	
	# print starting point
	println("Starting priority list:")
	printPriorityList(priorityList)
	
	# calculate objective value for starting point
	objVal = objValLookup(priorityList)
	if objVal == nullObjVal
		objVal = simObjVal!(sim, priorityList)
	end
	println("starting objective value: ", objVal)
	bestObjVal = objVal # current best objective value
	bestPriorityList = copy(priorityList)
	
	# write starting point
	for (key, value) in logFileDict
		if key != "search"
			logFileDict[key] = "" # reset value
		end
	end
	logFileDict["iter"] = 0
	logFileDict["objVal"] = objVal
	logFileDict["bestObjVal"] = bestObjVal
	logFileDict["iterTimeSeconds"] = 0
	logFileWriteDlmLine!(logFile, logFileDict, priorityList)
	
	# local search - try re-inserting/swapping items in the
	# priority list into different positions
	println("starting local search")
	endPoint = nothing
	iter = 1 # count iterations performed
	startTime = time() # time each iteration
	while true
		for i = (numAmbs+numStations):-1:1, j = 1:numAmbs, move = ["insert", "swap"]
			# starting with i > numAmbs will try insert/swap new items into priority list first (instead of only moving existing items within list); should help with speed
			
			i != j || continue # insert / swap will make no difference, skip
			
			# change priority list, using re-insertion / swapping
			# if i <= numAmbs, move item in priorityList[i] to priorityList[j]
			# if i > numAmbs, try move slot from station i-numAmbs (if num spare slots > 0)
			print(move, ": i = ", i, ", j = ", j)
			priorityList = copy(bestPriorityList)
			if move == "insert"
				insert!(priorityList, i, j) # || continue
			elseif move == "swap"
				swap!(priorityList, i, j) # || continue
			end
			
			objVal = objValLookup(priorityList)
			usedLookup = false
			if objVal != nullObjVal
				# priority list already tried
				print("; done")
				usedLookup = true
			else
				objVal = simObjVal!(sim, priorityList)
			end
			merge!(logFileDict, simStats(sim)) # need to do this before resetting sim
			
			println("; objective value: ", objVal, " (best: ", bestObjVal, ")")
			
			usedMove = false
			if (sense == :max && bestObjVal < objVal) || (sense == :min && bestObjVal > objVal)
				bestObjVal = objVal
				bestPriorityList = priorityList
				usedMove = true
				
				println("made ", move,": i = ", i, ", j = ", j, "; new objective value: ", bestObjVal)
			end
			
			logFileDict["iter"] = iter
			logFileDict["move"] = move
			logFileDict["i"] = i
			logFileDict["j"] = j
			logFileDict["usedLookup"] = Int(usedLookup)
			logFileDict["usedMove"] = Int(usedMove)
			logFileDict["objVal"] = objVal
			logFileDict["bestObjVal"] = bestObjVal
			logFileDict["iterTimeSeconds"] = round(time() - startTime, digits = 2)
			logFileWriteDlmLine!(logFile, logFileDict, priorityList)
			
			iter += 1
			
			if usedMove || endPoint == nothing
				endPoint = (i,j,move) # set end point
			elseif endPoint == (i,j,move)
				println("Finished local search ", logFileDict["search"], " (of ", numSearches, "); solution:")
				printPriorityList(bestPriorityList)
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
	@assert(i != j)
	
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
	@assert(i != j)
	
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

# run
t = time()
repeatedLocalSearch()
println("Total runtime: ", round(time()-t, digits = 2), " seconds")
