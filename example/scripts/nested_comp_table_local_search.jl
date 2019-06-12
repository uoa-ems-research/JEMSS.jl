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
nestedCompTables = [] # leave empty (and set numSearches) if generating random nested compliance tables for random restarts
const numSearches = isempty(nestedCompTables) ? 1 : length(nestedCompTables) # number of local searches to perform
const nestedCompTableRng = Random.MersenneTwister(0) # useful for reproducing results, if using random restarts

# some parameter checks
@assert(isfile(configFilename))
@assert(isdir(outputFolder))
@assert(ObjVal <: Real)
@assert(typeof(nullObjVal) <: ObjVal)
@assert(sense == :min || sense == :max)
@assert((isempty(nestedCompTables) && numSearches >= 1) || (length(nestedCompTables) == numSearches))
@assert(all(nestedCompTable -> checkCompTableIsNested(nestedCompTable), nestedCompTables))

stationCapacities = [] # stationCapacities[i] gives ambulance holding capacity of station i; will populate after sim is initialised

# keep track of nested compliance tables tried, and their objective values
nestedCompTablesObjVal = Dict{NestedCompTable, ObjVal}()
function objValLookup(nestedCompTable::NestedCompTable)::ObjVal
	return get(nestedCompTablesObjVal, nestedCompTable, nullObjVal) # return objective value if found, otherwise nullObjVal
end

# for a completed simulation, calculate and return the objective function value
function objFn(sim::Simulation)::ObjVal
	@assert(sim.complete)
	return countCallsReachedInTime(sim) # countCallsReachedInTime is from JEMSS
end

# get objective value for sim, applying nestedCompTable
# only use this if nestedCompTablesObjVal does not already contain nestedCompTable
# mutates: sim, nestedCompTablesObjVal
function simObjVal!(sim::Simulation, nestedCompTable::NestedCompTable)::ObjVal
	@assert(!haskey(nestedCompTablesObjVal, nestedCompTable))
	reset!(sim)
	for i = 1:sim.numAmbs
		setAmbStation!(sim.ambulances[i], sim.stations[nestedCompTable[i]])
	end
	initCompTable!(sim, nestedCompTable)
	simulateToEnd!(sim)
	objVal = objFn(sim)
	nestedCompTablesObjVal[nestedCompTable] = objVal
	# do not reset sim before returning, as other stats may needed after calling this function
	return objVal
end

function printNestedCompTable(nestedCompTable::NestedCompTable)
	for i = 1:length(nestedCompTable)
		println("item ", i, ": station ", nestedCompTable[i])
	end
end

function simStats(sim::Simulation)
	# some sim statistics to write to log file
	stats = Dict{String,Any}()
	stats["totalAmbTravelTime"] = ""
	stats["totalAmbBusyTime"] = ""
	stats["avgResponseDurationMinutes"] = ""
	stats["callsReachedInTime"] = ""
	if sim.complete
		stats["totalAmbTravelTime"] = sum(amb -> amb.totalTravelTime, sim.ambulances)
		stats["totalAmbBusyTime"] = sum(amb -> amb.totalBusyTime, sim.ambulances)
		@assert(all(call -> call.responseDuration != nullTime, sim.calls))
		stats["avgResponseDurationMinutes"] = mean(call -> call.responseDuration, sim.calls) * 24 * 60
		stats["callsReachedInTime"] = sum(call -> call.responseDuration <= sim.targetResponseDurations[Int(call.priority)], sim.calls)
	end
	return stats
end

logFileHeader = ["search", "iter", "move", "i", "j", "usedLookup", "usedMove", "objVal", "bestObjVal", "totalAmbTravelTime", "totalAmbBusyTime", "avgResponseDurationMinutes", "callsReachedInTime", "iterTimeSeconds"] # ... and nested comp table
logFileDict = Dict{String,Any}([s => "" for s in logFileHeader])

function logFileWriteDlmLine!(logFile::IOStream, data::Dict{String,Any}, nestedCompTable::NestedCompTable)
	line = []
	for header in logFileHeader
		push!(line, data[header])
	end
	line = vcat(line, nestedCompTable)
	writeDlmLine!(logFile, line...)
	flush(logFile)
end

nestedCompTableSols = Vector{NestedCompTable}() # solution of each complete local search

# perform local search, starting at each of the nested compliance tables provided
function repeatedLocalSearch()
	println("Initialising simulation from config: ", configFilename)
	sim = initSim(configFilename; createBackup = false, doPrint = false)
	
	# change move up module index to comp_table
	sim.moveUpData.moveUpModule = compTableModule
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
	
	global nestedCompTables, nestedCompTableSols
	if isempty(nestedCompTables)
		println("generating $numSearches random nested compliance table(s)")
		for i = 1:numSearches
			nestedCompTable = makeRandNestedCompTable(sim.numAmbs, sim.numStations;
				stationCapacities = stationCapacities, rng = nestedCompTableRng)
			push!(nestedCompTables, nestedCompTable)
		end
	end
	
	for i = 1:numSearches
		logFileDict["search"] = i
		
		# perform local search until no improvements can be made
		println()
		println("Search ", i, " (of ", numSearches, ")")
		nestedCompTable = localSearch!(sim, nestedCompTables[i], logFile)
		
		# save best ambulance to station allocation found for this search iteration
		push!(nestedCompTableSols, nestedCompTable)
		
		# write solution to file
		writeDlmLine!(solFile, i, objValLookup(nestedCompTable), nestedCompTable...)
		flush(solFile)
	end
	
	# print out all results from each finished local search
	println()
	println("Nested compliance tables from completed local searches:")
	for i = 1:numSearches
		println()
		println("Iteration: ", i)
		nestedCompTable = nestedCompTableSols[i]
		printNestedCompTable(nestedCompTable)
		println("objective value = ", objValLookup(nestedCompTable))
	end
	
	close(solFile)
	close(logFile)
end

# Local search of nested compliance table.
# Move items in nested compliance table, making changes that optimise
# the objective function (objFn) for the given sense (:min or :max).
# Continue search until no improvement can be made.
# Mutates: sim, logFile
function localSearch!(sim::Simulation, nestedCompTable::NestedCompTable, logFile::IOStream)
	
	# shorthand
	numAmbs = sim.numAmbs
	numStations = sim.numStations
	
	# print starting point
	println("Starting nested comp table:")
	printNestedCompTable(nestedCompTable)
	
	# calculate objective value for starting point
	objVal = objValLookup(nestedCompTable)
	if objVal == nullObjVal
		objVal = simObjVal!(sim, nestedCompTable)
	end
	println("starting objective value: ", objVal)
	bestObjVal = objVal # current best objective value
	bestNestedCompTable = copy(nestedCompTable)
	
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
	logFileWriteDlmLine!(logFile, logFileDict, nestedCompTable)
	
	# local search - try re-inserting/swapping items in the
	# nested comp table into different positions
	println("starting local search")
	endPoint = nothing
	iter = 1 # count iterations performed
	startTime = time() # time each iteration
	while true
		for i = (numAmbs+numStations):-1:1, j = 1:numAmbs, move = ["insert", "swap"]
			# starting with i > numAmbs will try insert/swap new items into nested comp table first (instead of only moving existing items within table); should help with speed
			
			i != j || continue # insert / swap will make no difference, skip
			
			# change nested comp table, using re-insertion / swapping
			# if i <= numAmbs, move item in nestedCompTable[i] to nestedCompTable[j]
			# if i > numAmbs, try move slot from station i-numAmbs (if num spare slots > 0)
			print(move, ": i = ", i, ", j = ", j)
			nestedCompTable = copy(bestNestedCompTable)
			if move == "insert"
				insert!(nestedCompTable, i, j) # || continue
			elseif move == "swap"
				swap!(nestedCompTable, i, j) # || continue
			end
			
			objVal = objValLookup(nestedCompTable)
			usedLookup = false
			if objVal != nullObjVal
				# nested comp table already tried
				print("; done")
				usedLookup = true
			else
				objVal = simObjVal!(sim, nestedCompTable)
			end
			merge!(logFileDict, simStats(sim)) # need to do this before resetting sim
			
			println("; objective value: ", objVal, " (best: ", bestObjVal, ")")
			
			usedMove = false
			if (sense == :max && bestObjVal < objVal) || (sense == :min && bestObjVal > objVal)
				bestObjVal = objVal
				bestNestedCompTable = nestedCompTable
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
			logFileWriteDlmLine!(logFile, logFileDict, nestedCompTable)
			
			iter += 1
			
			if usedMove || endPoint == nothing
				endPoint = (i,j,move) # set end point
			elseif endPoint == (i,j,move)
				println("Finished local search ", logFileDict["search"], " (of ", numSearches, "); solution:")
				printNestedCompTable(bestNestedCompTable)
				return bestNestedCompTable
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

# Insert item at index i into index j in nested comp table.
# If i > numAmbs, insertion is done with station index i-numAmbs.
# Return true if insertion completed, false otherwise.
function insert!(nestedCompTable::NestedCompTable, i::Int, j::Int)
	@assert(i != j)
	
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

# Swap items at index i and index j in nested comp table.
# If i > numAmbs, swap is done with station index i-numAmbs.
# Return true if swap completed, false otherwise.
function swap!(nestedCompTable::NestedCompTable, i::Int, j::Int)
	@assert(i != j)
	
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

# run
t = time()
repeatedLocalSearch()
println("Total runtime: ", round(time()-t, digits = 2), " seconds")
