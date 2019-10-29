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
using StatsFuns
using Base.Iterators: Stateful

const StationsNumAmbs = Vector{Int} # type alias

function deploymentLocalSearch!(sim::Simulation, deployments::Vector{Deployment}; outputFolder::String = "")
	isdir(outputFolder) || mkdir(outputFolder)
	@assert(isdir(outputFolder))
	@assert(!isempty(deployments))
	
	# parameters
	global solFilename = "$outputFolder/solutions.csv" # final solutions from each local search
	global deploymentsOutputFilename = "$outputFolder/deployments.csv" # solutions, stored as deployments
	global logFilename = "$outputFolder/log.csv" # log search progress to file
	global sense = :max # :min or :max; direction of optimisation for objective function
	global conf = 0.95 # statistical confidence level
	global doPrint = true
	
	# some parameter checks
	@assert(sense == :min || sense == :max)
	@assert(0 <= conf < 1) # should be between 0 and 1, but not equal to 1 otherwise confidence interval is infinite
	
	global nullObjVal = MeanAndHalfWidth(NaN,NaN)
	
	repeatedLocalSearch!(sim, deployments)
end

# keep track of station ambulance counts tried, and their objective values
global stationsNumAmbsObjVal = Dict{StationsNumAmbs, MeanAndHalfWidth}()
global stationsNumAmbsStats = Dict{StationsNumAmbs, Dict{String, Any}}()
function objValLookup(stationsNumAmbs::StationsNumAmbs)::MeanAndHalfWidth
	return get(stationsNumAmbsObjVal, stationsNumAmbs, nullObjVal) # return objective value if found, otherwise nullObjVal
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

# apply stationsNumAmbs to sim replications, and return objective value (and half-width) and whether lookup was used
# uses look-up if possible
# mutates: sim, stationsNumAmbsObjVal
function simObjVal!(sim::Simulation, stationsNumAmbs::StationsNumAmbs)::Tuple{MeanAndHalfWidth,Bool}
	global stationsNumAmbsObjVal
	objVal = objValLookup(stationsNumAmbs)
	if objVal != nullObjVal return (objVal, true) end
	for rep in sim.reps
		reset!(sim)
		applyStationsNumAmbs!(sim, stationsNumAmbs)
		simulateRep!(sim, rep)
		@assert(getStationsNumAmbs(sim) == stationsNumAmbs) # check that sim had the current deployment applied
	end
	objVal = objFn(sim)
	stationsNumAmbsObjVal[stationsNumAmbs] = objVal
	return (objVal, false)
end

# save the locally optimal station ambulance counts found for each search
global stationsNumAmbsSols = Vector{StationsNumAmbs}() # stationsNumAmbsSols[i] is solution for deployments[i]

# print number of ambulances at each station
function printStationsNumAmbs(stationsNumAmbs::StationsNumAmbs)
	for i = 1:length(stationsNumAmbs)
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

global solFileHeader = ["search", "objVal", "objValHalfWidth"] # ... and stationsNumAmbs
global logFileHeader = vcat(["search", "iter", "i", "j", "usedLookup", "usedMove", "searchDurationSeconds", "objVal", "objValHalfWidth", "bestObjVal"], sort(collect(keys(simStatsEmpty)))) # ... and stationsNumAmbs
global logFileDict = Dict{String,Any}([s => "" for s in logFileHeader])

# mutates: file
function fileWriteDlmLine!(file::IOStream, fileHeader::Vector{String}, data::Dict{String,T}, stationsNumAmbs::StationsNumAmbs) where T <: Any
	line = []
	for header in fileHeader
		push!(line, data[header])
	end
	line = vcat(line, stationsNumAmbs)
	writeDlmLine!(file, line...)
	flush(file)
end

# perform local search, starting at each of the deployments provided
# mutates: sim
function repeatedLocalSearch!(sim::Simulation, deployments::Vector{Deployment})
	global doPrint
	
	# open files for writing solution
	global solFilename, logFilename
	solFile = open(solFilename, "w")
	logFile = open(logFilename, "w")
	
	# write misc data to files
	global warmUpDuration, periodDuration, conf
	numSearches = length(deployments)
	miscTable = Table("miscData", ["numAmbs", "numStations", "numSearches", "conf"];
		rows = [[sim.numAmbs, sim.numStations, numSearches, conf]]) # note that numCalls may vary between sim replications
	repsTable = Table("repsData", ["numReps", "warmUpDuration", "periodDuration"];
		rows = [[sim.numReps, warmUpDuration, periodDuration]])
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
		
		# perform local search until no improvements can be made
		doPrint && println()
		doPrint && println("Search ", i, " (of ", numSearches, ")")
		doPrint && printStationsNumAmbs(stationsNumAmbs)
		stationsNumAmbs = localSearch!(sim, stationsNumAmbs, logFile)
		doPrint && println("Finished local search $i (of $numSearches); solution:")
		doPrint && printStationsNumAmbs(stationsNumAmbs)
		
		# save best station ambulance count found for this search
		global stationsNumAmbsSols
		push!(stationsNumAmbsSols, stationsNumAmbs)
		
		# write solution to file
		objValMeanAndHalfWidth = objValLookup(stationsNumAmbs)
		solFileDict = Dict("search" => i, "objVal" => objValMeanAndHalfWidth.mean, "objValHalfWidth" => objValMeanAndHalfWidth.halfWidth)
		fileWriteDlmLine!(solFile, solFileHeader, solFileDict, stationsNumAmbs)
		
		global deploymentsOutputFilename
		writeDeploymentsFile(deploymentsOutputFilename, map(stationsNumAmbsToDeployment, stationsNumAmbsSols), sim.numStations)
	end
	
	# print out all results from each finished local search
	doPrint && println()
	doPrint && println("Station ambulance counts from completed local searches:")
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

# Local search of ambulance deployment, starting at stationsNumAmbs.
# Move ambulances between stations, making changes that optimise
# the objective function (objFn) for the given sense (:min or :max).
# Continue search until no improvement can be made.
# mutates: sim, logFile
function localSearch!(sim::Simulation, stationsNumAmbs::StationsNumAmbs, logFile::IOStream)::StationsNumAmbs
	
	global logFileHeader, logFileDict, stationsNumAmbsStats, sense, doPrint
	
	# track iterations
	iter = 1 # count iterations performed
	startTime = time()
	getSearchDuration() = round(time() - startTime, digits = 2)
	
	# calculate objective value for starting point
	objValMeanAndHalfWidth, usedLookup = simObjVal!(sim, stationsNumAmbs)
	(objVal, objValHalfWidth) = (objValMeanAndHalfWidth.mean, objValMeanAndHalfWidth.halfWidth)
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
	if !usedLookup stationsNumAmbsStats[stationsNumAmbs] = stats end
	merge!(logFileDict, stats)
	merge!(logFileDict, Dict("iter" => iter, "usedLookup" => Int(usedLookup), "searchDurationSeconds" => getSearchDuration(),
		"objVal" => objVal, "objValHalfWidth" => objValHalfWidth, "bestObjVal" => bestObjVal))
	fileWriteDlmLine!(logFile, logFileHeader, logFileDict, stationsNumAmbs)
	
	# local search
	doPrint && println("starting local search")
	endPoint = nothing
	iter += 1
	numStations = sim.numStations # shorthand
	while true
		for i = 1:numStations, j = [i+1:numStations; 1:i-1;]
			newBestFound = false
			if bestStationsNumAmbs[i] > 0
				# test moving an ambulance from station i to j
				
				doPrint && print("move: i = ", i, ", j = ", j)
				
				# change station ambulance counts, later keep if improvement made
				stationsNumAmbs = copy(bestStationsNumAmbs)
				stationsNumAmbs[i] -= 1
				stationsNumAmbs[j] += 1
				
				objValMeanAndHalfWidth, usedLookup = simObjVal!(sim, stationsNumAmbs)
				(objVal, objValHalfWidth) = (objValMeanAndHalfWidth.mean, objValMeanAndHalfWidth.halfWidth)
				if !usedLookup
					stats = stationsNumAmbsStats[stationsNumAmbs] = getSimStats(sim)
					merge!(logFileDict, stats)
				else
					merge!(logFileDict, stationsNumAmbsStats[stationsNumAmbs])
					doPrint && print("; done")
				end
				
				doPrint && println("; objective value: ", objVal, " (best: ", bestObjVal, ")")
				
				usedMove = false
				if (sense == :max && bestObjVal < objVal) || (sense == :min && bestObjVal > objVal)
					bestObjVal = objVal
					# keep change (move ambulance from station i to j)
					bestStationsNumAmbs = stationsNumAmbs
					usedMove = true
					newBestFound = true
					
					doPrint && println("moved: i = ", i, ", j = ", j, "; new objective value: ", bestObjVal)
				end
				
				merge!(logFileDict, Dict("iter" => iter, "i" => i, "j" => j, "usedLookup" => Int(usedLookup), "usedMove" => Int(usedMove),
					"searchDurationSeconds" => getSearchDuration(), "objVal" => objVal, "objValHalfWidth" => objValHalfWidth, "bestObjVal" => bestObjVal))
				fileWriteDlmLine!(logFile, logFileHeader, logFileDict, stationsNumAmbs)
				
				iter += 1
			end
			
			if newBestFound || endPoint == nothing
				endPoint = (i,j)
			elseif endPoint == (i,j)
				doPrint && println("Found local optimum.")
				return bestStationsNumAmbs
			end
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
sim = initSim(configFilename, createBackup = true)
@assert(sim.numReps >= 2, "Need at least two simulation replications (sim.reps) for statistics.")

# set sim stats capturing
setSimStatsCapture!(sim, Stateful(periodDuration), warmUpDuration)
periodEndTime = sim.startTime + warmUpDuration + periodDuration
for rep in sim.reps
	# statistics should stop collecting before last call arrives (if not earlier), to avoid cool-down period
	@assert(periodEndTime <= rep.calls[end].arrivalTime)
end

# generate deployments
deployments = makeRandDeployments(sim, numSearches; rng = deploymentRng)

# run
deploymentLocalSearch!(sim, deployments, outputFolder = outputFolder)
println("total runtime: ", round(time()-t, digits = 2), " seconds")
