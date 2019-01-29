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

# Run a local search algorithm to optimise ambulance deployment, given some objective function (objFn).
# The local search is simple hill climbing, and the neighbourhood of a solution
# (deployment) involves moving one ambulance from its station to another station.
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
deployments = [] # leave empty (and set numSearches) if generating random deployments for random restarts
const numSearches = isempty(deployments) ? 1 : length(deployments) # number of local searches to perform
deploymentRng = Random.MersenneTwister(0) # useful for reproducing results, if using random restarts

# some parameter checks
@assert(isfile(configFilename))
@assert(isdir(outputFolder))
@assert(ObjVal <: Real)
@assert(typeof(nullObjVal) <: ObjVal)
@assert(sense == :min || sense == :max)
@assert((isempty(deployments) && numSearches >= 1) || (length(deployments) == numSearches))

const StationsNumAmbs = Vector{Int} # type alias

# keep track of station ambulance counts tried, and their objective values
stationsNumAmbsObjVal = Dict{StationsNumAmbs, ObjVal}()
function objValLookup(stationsNumAmbs::StationsNumAmbs)::ObjVal
	return get(stationsNumAmbsObjVal, stationsNumAmbs, nullObjVal) # return objective value if found, otherwise nullObjVal
end

# for a completed simulation, calculate and return the objective function value
function objFn(sim::Simulation)::ObjVal
	@assert(sim.complete)
	return countCallsReachedInTime(sim) # countCallsReachedInTime is from JEMSS
end

# get objective value for sim, applying stationsNumAmbs
# only use this if stationsNumAmbsObjVal does not already contain stationsNumAmbs
# mutates: sim, stationsNumAmbsObjVal
function simObjVal!(sim::Simulation, stationsNumAmbs::StationsNumAmbs)::ObjVal
	@assert(!haskey(stationsNumAmbsObjVal, stationsNumAmbs))
	resetSim!(sim)
	applyStationsNumAmbs!(sim, stationsNumAmbs)
	simulateToEnd!(sim)
	objVal = objFn(sim)
	stationsNumAmbsObjVal[stationsNumAmbs] = objVal
	# do not reset sim before returning, as other stats may needed after calling this function
	return objVal
end

# save the locally optimal station ambulance counts found for each search
stationsNumAmbsSols = Vector{StationsNumAmbs}() # stationsNumAmbsSols[i] is solution for deployments[i]

# print number of ambulances at each station
function printStationsNumAmbs(stationsNumAmbs::StationsNumAmbs)
	for i = 1:length(stationsNumAmbs)
		println("station ", i, ": ", stationsNumAmbs[i], " ambulances")
	end
end

function simStats(sim::Simulation)
	# some sim statistics to write to log file
	stats = Dict{String,Any}()
	stats["totalAmbTravelTime"] = ""
	stats["totalAmbBusyTime"] = ""
	stats["avgResponseTimeMinutes"] = ""
	stats["callsReachedInTime"] = ""
	if sim.complete
		stats["totalAmbTravelTime"] = sum(amb -> amb.totalTravelTime, sim.ambulances)
		stats["totalAmbBusyTime"] = sum(amb -> amb.totalBusyTime, sim.ambulances)
		@assert(all(call -> call.responseTime != nullTime, sim.calls))
		stats["avgResponseTimeMinutes"] = mean(call -> call.responseTime, sim.calls) * 24 * 60
		stats["callsReachedInTime"] = sum(call -> call.responseTime <= sim.targetResponseTimes[Int(call.priority)], sim.calls)
	end
	return stats
end

logFileHeader = ["search", "iter", "i", "j", "usedLookup", "usedMove", "objVal", "bestObjVal", "totalAmbTravelTime", "totalAmbBusyTime", "avgResponseTimeMinutes", "callsReachedInTime", "iterTimeSeconds"] # ... and stationsNumAmbs
logFileDict = Dict{String,Any}([s => "" for s in logFileHeader])

function logFileWriteDlmLine!(logFile::IOStream, data::Dict{String,Any}, stationsNumAmbs::StationsNumAmbs)
	line = []
	for header in logFileHeader
		push!(line, data[header])
	end
	line = vcat(line, stationsNumAmbs)
	writeDlmLine!(logFile, line...)
	flush(logFile)
end

# perform local search, starting at each of the deployments provided
function repeatedLocalSearch()
	println("Initialising simulation from config: ", configFilename)
	sim = initSim(configFilename, createBackup = true, doPrint = false)
	
	# open files for writing solution
	solFile = open(solFilename, "w")
	logFile = open(logFilename, "w")
	
	# write file headers
	writeDlmLine!(solFile, "search", "objVal", "solution:", [1:sim.numStations;]...)
	writeDlmLine!(logFile, logFileHeader..., ["station_$i numAmbs" for i = 1:sim.numStations]...)
	flush(solFile)
	flush(logFile)
	
	global deployments
	if isempty(deployments)
		println("generating $numSearches random deployment(s)")
		deployments = makeRandDeployments(sim, numSearches; rng = deploymentRng)
	end
	
	for i = 1:numSearches
		
		logFileDict["search"] = i
		
		# assign ambulances to stations
		applyDeployment!(sim, deployments[i])
		
		# perform local search until no improvements can be made
		println()
		println("Search ", i, " (of ", numSearches, ")")
		stationsNumAmbs = localSearch!(sim, stationsNumAmbsObjVal, logFile)
		
		# save best station ambulance count found for this search iteration
		push!(stationsNumAmbsSols, stationsNumAmbs)
		
		# write solution to file
		writeDlmLine!(solFile, i, objValLookup(stationsNumAmbs), "", stationsNumAmbs...)
		flush(solFile)
		
		resetSim!(sim)
	end
	
	# print out all results from each finished local search
	println()
	println("Station ambulance counts from completed local searches:")
	for i = 1:numSearches
		println()
		println("Iteration: ", i)
		stationsNumAmbs = stationsNumAmbsSols[i]
		printStationsNumAmbs(stationsNumAmbs)
		println("objective value = ", objValLookup(stationsNumAmbs))
	end
	
	close(solFile)
	close(logFile)
end

# Local search of ambulance deployment.
# Move ambulances between stations, making changes that optimise
# the objective function (objFn) for the given sense (:min or :max).
# Continue search until no improvement can be made.
# mutates: sim, stationsNumAmbsObjVal, logFile
function localSearch!(sim::Simulation, stationsNumAmbsObjVal::Dict{StationsNumAmbs,ObjVal}, logFile::IOStream)::StationsNumAmbs
	
	# shorthand
	numStations = sim.numStations
	
	# count number of ambs at each station
	stationsNumAmbs = getStationsNumAmbs(sim)
	
	# print number of ambulances at each station
	printStationsNumAmbs(stationsNumAmbs)
	
	# calculate objective value for starting point
	objVal = objValLookup(stationsNumAmbs)
	if objVal == nullObjVal
		objVal = simObjVal!(sim, stationsNumAmbs)
	end
	println("starting objective value: ", objVal)
	bestObjVal = objVal # current best objective value
	bestStationsNumAmbs = copy(stationsNumAmbs)
	
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
	logFileWriteDlmLine!(logFile, logFileDict, stationsNumAmbs)
	
	# local search
	println("starting local search")
	endPoint = nothing
	iter = 1 # count iterations performed
	startTime = time() # time each iteration
	while true
		for i = 1:numStations, j = [i+1:numStations; 1:i-1;]
			newBestFound = false
			if bestStationsNumAmbs[i] > 0
				# test moving an ambulance from station i to j
				
				print("move: i = ", i, ", j = ", j)
				
				# resetSim!(sim) # only resets if sim was used
				
				# change station ambulance counts, later keep if improvement made
				stationsNumAmbs = copy(bestStationsNumAmbs)
				stationsNumAmbs[i] -= 1
				stationsNumAmbs[j] += 1
				
				objVal = objValLookup(stationsNumAmbs)
				usedLookup = false
				if objVal == nullObjVal
					objVal = simObjVal!(sim, stationsNumAmbs)
				else
					print("; done")
					usedLookup = true
				end
				merge!(logFileDict, simStats(sim)) # need to do this before resetting sim
				
				println("; objective value: ", objVal, " (best: ", bestObjVal, ")")
				
				usedMove = false
				if (sense == :max && bestObjVal < objVal) || (sense == :min && bestObjVal > objVal)
					bestObjVal = objVal
					# keep change (move ambulance from station i to j)
					bestStationsNumAmbs = stationsNumAmbs
					usedMove = true
					newBestFound = true
					
					println("moved: i = ", i, ", j = ", j, "; new objective value: ", bestObjVal)
				end
				
				logFileDict["iter"] = iter
				logFileDict["i"] = i
				logFileDict["j"] = j
				logFileDict["usedLookup"] = Int(usedLookup)
				logFileDict["usedMove"] = Int(usedMove)
				logFileDict["objVal"] = objVal
				logFileDict["bestObjVal"] = bestObjVal
				logFileDict["iterTimeSeconds"] = round(time() - startTime, digits = 2)
				logFileWriteDlmLine!(logFile, logFileDict, stationsNumAmbs)
				
				iter += 1
			end
			
			if newBestFound || endPoint == nothing
				endPoint = (i,j)
			elseif endPoint == (i,j)
				println("Finished local search ", logFileDict["search"], " (of ", numSearches, "); solution:")
				printStationsNumAmbs(bestStationsNumAmbs)
				return bestStationsNumAmbs
			end
		end
	end
end

# run
t = time()
repeatedLocalSearch()
println("total runtime: ", round(time()-t, digits = 2), " seconds")
