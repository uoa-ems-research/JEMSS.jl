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

# write sim object and output files

function writeAmbsFile(filename::String, ambulances::Vector{Ambulance}; writeOutputFields::Bool = false)
	header = ["index", "stationIndex", "class"]
	row1(a::Ambulance) = [a.index, a.stationIndex, Int(a.class)]
	
	if writeOutputFields
		counts = (:numCallsTreated, :numCallsTransported, :numDispatches, :numDispatchesFromStation, :numDispatchesOnRoad, :numDispatchesOnFree, :numRedispatches, :numMoveUps, :numMoveUpsFromStation, :numMoveUpsOnRoad, :numMoveUpsOnFree, :numMoveUpsReturnToPrevStation)
		countHeaders = [string(c) for c in counts]
		
		statuses = (setdiff(instances(AmbStatus), (ambNullStatus,))..., instances(AmbStatusSet)...)
		travelStatuses = (ambStatusSets[ambTravelling]..., instances(AmbStatusSet)...)
		statusDurationHeaders = [string("statusDurations_", string(s)) for s in statuses]
		statusDistanceHeaders = [string("statusDistances_", string(s)) for s in travelStatuses]
		
		header = vcat(header, countHeaders, statusDurationHeaders, statusDistanceHeaders)
		row2(a::Ambulance) = vcat([getfield(a,c) for c in counts], [a.statusDurations[s] for s in statuses], [a.statusDistances[s] for s in travelStatuses])
		
		# skipped: statusTransitionCounts (array)
	end
	
	row(a::Ambulance) = writeOutputFields ? vcat(row1(a), row2(a)) : row1(a)
	table = Table("ambulances", header; rows = [row(a) for a in ambulances])
	writeTablesToFile(filename, table)
end

function writeArcsFile(filename::String, arcs::Vector{Arc}, travelTimes::Array{Float,2}, arcForm::String)
	@assert(arcForm == "directed" || arcForm == "undirected")
	numModes = size(travelTimes,1)
	miscTable = Table("miscData", ["arcForm", "numModes"]; rows = [[arcForm, numModes]])
	mainHeaders = ["index", "fromNode", "toNode", "distance", ["mode_$i" for i = 1:numModes]...]
	fieldNames = setdiff(collect(keys(arcs[1].fields)), mainHeaders) # assume fields are same for all arcs
	arcsTable = Table("arcs", [mainHeaders..., fieldNames...];
		rows = [vcat(a.index, a.fromNodeIndex, a.toNodeIndex, a.distance, travelTimes[:,a.index]..., [a.fields[f] for f in fieldNames]...) for a in arcs])
	writeTablesToFile(filename, [miscTable, arcsTable])
end

function writeCallsFile(filename::String, startTime::Float, calls::Vector{Call}; writeOutputFields::Bool = false)
	@assert(length(calls) >= 1)
	miscTable = Table("miscData", ["startTime"]; rows = [[startTime]])
	
	header = ["index", "priority", "x", "y", "arrivalTime", "dispatchDelay", "onSceneDuration", "transport", "hospitalIndex", "handoverDuration"]
	row1(c::Call) = [c.index, Int(c.priority), c.location.x, c.location.y, c.arrivalTime, c.dispatchDelay, c.onSceneDuration, Int(c.transport), c.hospitalIndex, c.handoverDuration]
	
	if writeOutputFields
		header = vcat(header, ["dispatchTime", "ambArrivalTime", "hospitalArrivalTime", "numBumps", "wasQueued", "ambDispatchLoc.x", "ambDispatchLoc.y", "ambStatusBeforeDispatch", "chosenHospitalIndex", "queuedDuration", "bumpedDuration", "waitingForAmbDuration", "responseDuration", "ambGoingToCallDuration", "transportDuration", "serviceDuration"])
		row2(c::Call) = [c.dispatchTime, c.ambArrivalTime, c.hospitalArrivalTime, c.numBumps, Int(c.wasQueued), c.ambDispatchLoc.x, c.ambDispatchLoc.y, string(c.ambStatusBeforeDispatch), c.chosenHospitalIndex, c.queuedDuration, c.bumpedDuration, c.waitingForAmbDuration, c.responseDuration, c.ambGoingToCallDuration, c.transportDuration, c.serviceDuration]
	end
	
	row(c::Call) = writeOutputFields ? vcat(row1(c), row2(c)) : row1(c)
	callsTable = Table("calls", header; rows = [row(c) for c in calls])
	writeTablesToFile(filename, [miscTable, callsTable])
end

function writeDemandFile(filename::String, demand::Demand)
	demandRastersTable = Table("demandRasters", ["rasterIndex", "rasterFilename"];
		rows = [[i, demand.rasterFilenames[i]] for i = 1:demand.numRasters])
	
	demandModesTable = Table("demandModes", ["modeIndex", "rasterIndex", "priority", "arrivalRate"];
		rows = [[m.index, m.rasterIndex, string(m.priority), m.arrivalRate] for m in demand.modes])
	@assert(all(i -> demand.modes[i].index == i, 1:demand.numModes))
	
	# demand sets table
	dml = demand.modeLookup # shorthand
	@assert(size(dml) == (demand.numSets, numPriorities)) # should have value for each combination of demand mode and priority
	demandSetsTable = Table("demandSets", ["setIndex", "modeIndices"];
		rows = [[i, hcat(dml[i,:]...)] for i = 1:demand.numSets])
	
	# demand sets timing table
	startTimes = demand.setsStartTimes # shorthand
	setIndices = demand.setsTimeOrder # shorthand
	demandSetsTimingTable = Table("demandSetsTiming", ["startTime", "setIndex"];
		rows = [[startTimes[i], setIndices[i]] for i = 1:length(startTimes)])
	
	writeTablesToFile(filename, [demandRastersTable, demandModesTable, demandSetsTable, demandSetsTimingTable])
end

function writeDemandCoverageFile(filename::String, demandCoverage::DemandCoverage)
	dc = demandCoverage # shorthand
	coverTimesTable = Table("coverTimes", ["demandPriority", "coverTime"];
		rows = [[string(priority), coverTime] for (priority, coverTime) in dc.coverTimes])
	
	demandRasterCellNumPointsTable = Table("demandRasterCellNumPoints", ["rows", "cols"];
		rows = [[dc.rasterCellNumRows, dc.rasterCellNumCols]])
	
	writeTablesToFile(filename, [coverTimesTable, demandRasterCellNumPointsTable])
end

function writeHospitalsFile(filename::String, hospitals::Vector{Hospital}; writeOutputFields::Bool = false)
	header = ["index", "x", "y"]
	row1(h::Hospital) = [h.index, h.location.x, h.location.y]
	
	if writeOutputFields
		header = vcat(header, ["numCalls"])
		row2(h::Hospital) = [h.numCalls]
	end
	
	row(h::Hospital) = writeOutputFields ? vcat(row1(h), row2(h)) : row1(h)
	table = Table("hospitals", header, rows = [row(h) for h in hospitals])
	writeTablesToFile(filename, table)
end

function writeMapFile(filename::String, map::Map)
	table = Table("map", ["xMin", "xMax", "yMin", "yMax", "xScale", "yScale"];
		rows = [[map.xMin, map.xMax, map.yMin, map.yMax, map.xScale, map.yScale]])
	writeTablesToFile(filename, table)
end

function writeNodesFile(filename::String, nodes::Vector{Node})
	mainHeaders = ["index", "x", "y", "offRoadAccess"]
	fieldNames = setdiff(collect(keys(nodes[1].fields)), mainHeaders) # assume fields are same for all nodes
	table = Table("nodes", [mainHeaders..., fieldNames...];
		rows = [[n.index, n.location.x, n.location.y, Int(n.offRoadAccess), [n.fields[f] for f in fieldNames]...] for n in nodes])
	writeTablesToFile(filename, table)
end

function writeRedispatchFile(filename::String, redispatch::Redispatch)
	miscTable = Table("miscData", ["allowRedispatch"], rows = [[Int(redispatch.allow)]])
	conditionsTable = Table("redispatchConditions", ["fromCallPriority", "toCallPriority", "allowRedispatch"],
		rows = [[[string(p1), string(p2), Int(redispatch.conditions[Int(p1),Int(p2)])] for p1 in priorities, p2 in priorities]...])
	writeTablesToFile(filename, [miscTable, conditionsTable])
end

function writeRNetTravelsFile(filename::String, rNetTravels::Vector{NetTravel})
	n = length(rNetTravels)
	@assert(all(i -> rNetTravels[i].isReduced, 1:n))
	@assert(all(i -> rNetTravels[i].modeIndex == i, 1:n))
	# save only some field values to file
	rNetTravelsSave = [NetTravel(true) for i = 1:n]
	for i = 1:n, fname in (:modeIndex, :arcTimes, :arcDists, :fadjList, :spFadjIndex, :spNodePairArcIndex, :spFadjArcList)
		setfield!(rNetTravelsSave[i], fname, getfield(rNetTravels[i], fname))
	end
	serializeToFile(filename, rNetTravelsSave)
end

function writePrioritiesFile(filename::String, targetResponseDurations::Vector{Float})
	table = Table("priorities", ["priority", "name", "targetResponseDuration"];
		rows = [[i, string(Priority(i)), targetResponseDurations[i]] for i = 1:length(targetResponseDurations)])
	writeTablesToFile(filename, table)
end

function writePriorityListFile(filename::String, priorityList::PriorityList)
	numAmbs = length(priorityList)
	table = Table("priorityList", ["item", "stationIndex"];
		cols = [collect(1:numAmbs), priorityList])
	writeTablesToFile(filename, table)
end

function writePriorityListsFile(filename::String, priorityLists::Vector{PriorityList})
	n = length(priorityLists)
	numAmbs = length(priorityLists[1])
	@assert(all(priorityList -> length(priorityList) == numAmbs, priorityLists))
	table = Table("priorityLists", ["item", ["priorityList_$i stationIndex" for i = 1:n]...];
		cols = [collect(1:numAmbs), priorityLists...])
	writeTablesToFile(filename, table)
end

function writeStationsFile(filename::String, stations::Vector{Station})
	table = Table("stations", ["index", "x", "y", "capacity"];
		rows = [[s.index, s.location.x, s.location.y, s.capacity] for s in stations])
	writeTablesToFile(filename, table)
end

function writeTravelFile(filename::String, travel::Travel)
	# travel modes table
	travelModesTable = Table("travelModes", ["travelModeIndex", "offRoadSpeed"];
		rows = [[tm.index, tm.offRoadSpeed] for tm in travel.modes])
	
	# travel sets table
	tml = travel.modeLookup # shorthand
	@assert(size(tml) == (length(travel.modes), numPriorities)) # should have value for each combination of travel mode and priority
	travelSetsTable = Table("travelSets", ["travelSetIndex", "priority", "travelModeIndex"];
		rows = [[[i, string(Priority(j)), tml[i,j]] for i = 1:size(tml,1), j = 1:size(tml,2)]...])
	
	# travel sets timing table
	startTimes = travel.setsStartTimes # shorthand
	setIndices = travel.setsTimeOrder # shorthand
	travelSetsTimingTable = Table("travelSetsTiming", ["startTime", "travelSetIndex"];
		rows = [[startTimes[i], setIndices[i]] for i = 1:length(startTimes)])
	
	writeTablesToFile(filename, [travelModesTable, travelSetsTable, travelSetsTimingTable])
end

function writeZhangIpParamsFile(filename::String, zhangIpData::ZhangIpData)
	# shorthand
	zid = zhangIpData
	@unpack marginalBenefits, stationCapacities = zid
	numStations = length(marginalBenefits)
	@assert(numStations == length(stationCapacities))
	
	miscTable = Table("miscParams", ["travelTimeCost", "onRoadMoveUpDiscountFactor", "regretTravelTimeThreshold", "expectedHospitalHandoverDuration"];
		rows = [[zid.travelTimeCost, zid.onRoadMoveUpDiscountFactor, zid.regretTravelTimeThreshold, zid.expectedHospitalHandoverDuration]])
	
	stationCapacitiesTable = Table("stationCapacities", ["stationIndex", "capacity"];
		rows = [[i, stationCapacities[i]] for i = 1:numStations])
	
	maxNumSlots = maximum(stationCapacities)
	mb = Array{Any,2}(undef, numStations, maxNumSlots)
	mb[:] .= ""
	for i = 1:numStations, j = 1:length(marginalBenefits[i])
		mb[i,j] = marginalBenefits[i][j]
	end
	marginalBenefitsTable = Table("stationMarginalBenefits", ["stationIndex", ["slot_$i" for i = 1:maxNumSlots]...];
		rows = [[i, mb[i,:]...] for i = 1:numStations])
	
	writeTablesToFile(filename, [miscTable, stationCapacitiesTable, marginalBenefitsTable])
end

# opens output files for writing during simulation
# note: should have field sim.resim.use set/fixed before calling this function
function openOutputFiles!(sim::Simulation)
	if !sim.writeOutput; return; end
	
	outputFilePath(name::String) = sim.outputFiles[name].path
	
	# create output path if it does not already exist
	if !isdir(sim.outputPath)
		mkdir(sim.outputPath)
	end
	
	# open output files for writing
	# currently, only need to open events file
	if !sim.resim.use # otherwise, existing events file is used for resimulation
		sim.outputFiles["events"].iostream = open(outputFilePath("events"), "w")
		sim.eventsFileIO = sim.outputFiles["events"].iostream # shorthand
		
		# save checksum of input files
		inputFiles = sort([name for (name, file) in sim.inputFiles])
		fileChecksumStrings = [string("'", sim.inputFiles[name].checksum, "'") for name in inputFiles]
		writeTablesToFile!(sim.eventsFileIO, Table("inputFiles", ["name", "checksum"]; cols = [inputFiles, fileChecksumStrings]))
		
		# save events with a key, to reduce file size
		eventForms = instances(EventForm)
		eventKeys = [Int(eventForm) for eventForm in eventForms]
		eventNames = [string(eventForm) for eventForm in eventForms]
		writeTablesToFile!(sim.eventsFileIO, Table("eventDict", ["key", "name"]; cols = [eventKeys, eventNames]))
		
		# write events table name and header
		writeDlmLine!(sim.eventsFileIO, "events")
		writeDlmLine!(sim.eventsFileIO, "index", "parentIndex", "time", "eventKey", "ambIndex", "callIndex", "stationIndex")
	end
end

function writeOutputFiles(sim::Simulation)
	writeMiscOutputFiles(sim)
	if !isempty(sim.reps)
		writeStatsFiles(sim, sim.reps)
	else
		writeStatsFiles(sim)
	end
end

function writeMiscOutputFiles(sim::Simulation)
	outputFileKeys = keys(sim.outputFiles)
	outputFilePath(name::String) = sim.outputFiles[name].path
	if in("ambulances", outputFileKeys) writeAmbsFile(outputFilePath("ambulances"), sim.ambulances; writeOutputFields = true) end
	if in("calls", outputFileKeys) writeCallsFile(outputFilePath("calls"), sim.startTime, sim.calls; writeOutputFields = true) end
	if in("hospitals", outputFileKeys) writeHospitalsFile(outputFilePath("hospitals"), sim.hospitals; writeOutputFields = true) end
end

function closeOutputFiles!(sim::Simulation)
	if !sim.writeOutput; return; end
	
	if !sim.resim.use
		writeDlmLine!(sim.eventsFileIO, "end")
		close(sim.eventsFileIO)
	end
end

function writeEventToFile!(sim::Simulation, event::Event)
	if !sim.writeOutput || sim.resim.use; return; end
	
	writeDlmLine!(sim.eventsFileIO, event.index, event.parentIndex, @sprintf("%0.5f", event.time), Int(event.form), event.ambIndex, event.callIndex, event.stationIndex)
	
	# flush(sim.eventsFileIO)
end

# write deployments to file
function writeDeploymentsFile(filename::String, deployments::Vector{Deployment}, numStations::Int)
	numAmbs = length(deployments[1])
	@assert(numStations >= maximum([maximum(d) for d in deployments]))
	numDeployments = length(deployments)
	
	miscTable = Table("miscData", ["numStations", "numDeployments"]; rows = [[numStations, numDeployments]])
	deploymentsTable = Table("deployments",
		["ambIndex", ["deployment_$i stationIndex" for i = 1:numDeployments]...];
		cols = [collect(1:numAmbs), deployments...])
	writeTablesToFile(filename, [miscTable, deploymentsTable])
end

# save batch mean response durations to file
function writeBatchMeanResponseDurationsFile(filename::String, batchMeanResponseDurations::Array{Float,2};
	batchTime = nullTime, startTime = nullTime, endTime = nullTime, responseDurationUnits::String = "minutes")
	@assert(batchTime != nullTime && startTime != nullTime && endTime != nullTime)
	x = batchMeanResponseDurations # shorthand
	(numRows, numCols) = size(x) # numRows = numSims, numCols = numBatches
	miscTable = Table("misc_data",
		["numSims", "numBatches", "batchTime", "startTime", "endTime", "response_duration_units"];
		rows=[[numRows, numCols, batchTime, startTime, endTime, responseDurationUnits]])
	avgBatchMeansTable = Table("avg_batch_mean_response_durations",
		["sim_index", "avg_batch_mean_response_duration", "standard_error"];
		rows = [[i, mean(x[i,:]), sem(x[i,:])] for i = 1:numRows])
	batchMeansTable = Table("batch_mean_response_durations",
		["batch_index", ["sim_$i" for i = 1:numRows]...];
		rows = [[i, x[:,i]...] for i = 1:numCols])
	writeTablesToFile(filename, [miscTable, avgBatchMeansTable, batchMeansTable])
end

# for batches in a single simulation replication
function writeStatsFiles(sim::Simulation)
	outputFileKeys = keys(sim.outputFiles)
	outputFilePath(name::String) = sim.outputFiles[name].path
	stats = sim.stats # shorthand
	if in("ambulancesStats", outputFileKeys) writeAmbsStatsFile(outputFilePath("ambulancesStats"), stats) end
	if in("callsStats", outputFileKeys) writeCallsStatsFile(outputFilePath("callsStats"), stats) end
	if in("hospitalsStats", outputFileKeys) writeHospitalsStatsFile(outputFilePath("hospitalsStats"), stats) end
	if in("stationsStats", outputFileKeys) writeStationsStatsFile(outputFilePath("stationsStats"), stats) end
	if in("statsDict", outputFileKeys) writeStatsDictFile(outputFilePath("statsDict"), getPeriodStatsList(sim)) end
end

# for independent simulation replications
function writeStatsFiles(sim::Simulation, reps::Vector{Simulation}; periodIndex::Int = nullIndex)
	outputFileKeys = keys(sim.outputFiles)
	outputFilePath(name::String) = sim.outputFiles[name].path
	periods = periodIndex == nullIndex ? getRepsPeriodStatsList(reps) : getRepsPeriodStatsList(reps, periodIndex = periodIndex)
	@assert(!isempty(periods))
	if in("ambulancesStats", outputFileKeys) writeAmbsStatsFile(outputFilePath("ambulancesStats"), periods) end
	if in("callsStats", outputFileKeys) writeCallsStatsFile(outputFilePath("callsStats"), periods) end
	if in("hospitalsStats", outputFileKeys) writeHospitalsStatsFile(outputFilePath("hospitalsStats"), periods) end
	if in("stationsStats", outputFileKeys) writeStationsStatsFile(outputFilePath("stationsStats"), periods) end
	if in("statsDict", outputFileKeys) writeStatsDictFile(outputFilePath("statsDict"), periods) end
end

function simStatsTimestampsTable(stats::SimStats)::Table
	return Table("timestamps", ["simStartTime", "warmUpEndTime", "lastCallArrivalTime", "simEndTime"];
		rows = [[stats.simStartTime, stats.warmUpEndTime, stats.lastCallArrivalTime, stats.simEndTime]])
end

function simStatsPeriodTable(period::SimPeriodStats)::Table
	return Table("period", ["startTime", "endTime", "duration"];
		rows = [[period.startTime, period.endTime, period.duration]])
end
function simStatsPeriodTable(periods::Vector{SimPeriodStats})::Table
	@assert(length(periods) >= 1)
	@unpack startTime, endTime, duration = periods[1]
	@assert(all(p -> isapprox(p.startTime, startTime), periods))
	@assert(all(p -> isapprox(p.endTime, endTime), periods))
	@assert(all(p -> isapprox(p.duration, duration), periods))
	return simStatsPeriodTable(periods[1])
end

function simStatsPeriodsTable(periods::Vector{SimPeriodStats})::Table
	periodsTable = Table("periods", ["periodIndex", "startTime", "endTime", "duration"];
		rows = [[i, p.startTime, p.endTime, p.duration] for (i,p) in enumerate(periods)])
	return periodsTable
end

function makeAmbsStatsTables(periods::Vector{SimPeriodStats})::Vector{Table}
	numAmbs = length(periods[1].ambulances) # shorthand
	@assert(all(p -> length(p.ambulances) == numAmbs, periods))
	
	ambulanceTables = Table[]
	counts = (:numCallsTreated, :numCallsTransported, :numDispatches, :numDispatchesFromStation, :numDispatchesOnRoad, :numDispatchesOnFree, :numRedispatches, :numMoveUps, :numMoveUpsFromStation, :numMoveUpsOnRoad, :numMoveUpsOnFree, :numMoveUpsReturnToPrevStation)
	statuses = (setdiff(instances(AmbStatus), (ambNullStatus,))..., instances(AmbStatusSet)...)
	travelStatuses = (sort(collect(ambStatusSets[ambTravelling]))..., instances(AmbStatusSet)...)
	countHeaders = [string(c) for c in counts]
	statusDurationHeaders = [string("duration_", string(s)) for s in statuses]
	statusDistanceHeaders = [string("distance_", string(s)) for s in travelStatuses]
	getAmb(period::SimPeriodStats, ambIndex::Int) = ambIndex == 0 ? period.ambulance : period.ambulances[ambIndex]
	header = vcat("periodIndex", countHeaders, statusDurationHeaders, statusDistanceHeaders)
	row(a::AmbulanceStats) = vcat([getfield(a,c) for c in counts], [a.statusDurations[s] for s in statuses], [a.statusDistances[s] for s in travelStatuses])
	# skipped: statusTransitionCounts
	for i = 0:numAmbs
		name = i == 0 ? "ambulance" : "ambulances[$i]"
		ambulanceTable = Table(name, header;
			rows = [vcat(j, row(getAmb(p,i))) for (j,p) in enumerate(periods)])
		push!(ambulanceTables, ambulanceTable)
	end
	
	return ambulanceTables
end

# for batches in a single simulation replication
function writeAmbsStatsFile(filename::String, stats::SimStats)
	# shorthand
	periods = stats.periods
	numAmbs = length(periods[1].ambulances)
	
	miscTable = Table("miscData", ["mode", "numAmbs"]; rows = [["batch", numAmbs]])
	timestampsTable = simStatsTimestampsTable(stats)
	periodsTable = simStatsPeriodsTable(periods)
	ambulanceTables = makeAmbsStatsTables(periods)
	writeTablesToFile(filename, [miscTable, timestampsTable, periodsTable, ambulanceTables...])
end

# for independent simulation replications
function writeAmbsStatsFile(filename::String, periods::Vector{SimPeriodStats})
	numAmbs = length(periods[1].ambulances) # shorthand
	@assert(all(p -> length(p.ambulances) == numAmbs, periods))
	
	miscTable = Table("miscData", ["mode", "numAmbs"]; rows = [["replication", numAmbs]])
	periodTable = simStatsPeriodTable(periods)
	ambulanceTables = makeAmbsStatsTables(periods)
	writeTablesToFile(filename, [miscTable, periodTable, ambulanceTables...])
end

function makeCallsStatsTables(periods::Vector{SimPeriodStats})::Vector{Table}
	fnames = setdiff(fieldnames(CallStats), (:callIndex,))
	callTable = Table("call", vcat("periodIndex", collect(string.(fnames)));
		rows = [vcat(i, [getfield(p.call, fname) for fname in fnames]) for (i,p) in enumerate(periods)])
	
	callPrioritiesTables = Table[]
	for priority in priorities
		table = Table("callPriorities[$priority]", vcat("periodIndex", collect(string.(fnames)));
			rows = [vcat(i, [getfield(p.callPriorities[priority], fname) for fname in fnames]) for (i,p) in enumerate(periods)])
		push!(callPrioritiesTables, table)
	end
	
	return [callTable, callPrioritiesTables...]
end

# for batches in a single simulation replication
function writeCallsStatsFile(filename::String, stats::SimStats)
	# shorthand
	periods = stats.periods
	numCalls = stats.captures[end].call.numCalls
	
	miscTable = Table("miscData", ["mode", "numCalls"]; rows = [["batch", numCalls]])
	timestampsTable = simStatsTimestampsTable(stats)
	periodsTable = simStatsPeriodsTable(periods)
	callsTables = makeCallsStatsTables(periods)
	writeTablesToFile(filename, [miscTable, timestampsTable, periodsTable, callsTables...])
end

# for independent simulation replications
function writeCallsStatsFile(filename::String, periods::Vector{SimPeriodStats})
	miscTable = Table("miscData", ["mode"]; rows = [["replication"]]) # note that number of calls may vary between replications
	periodTable = simStatsPeriodTable(periods)
	callsTables = makeCallsStatsTables(periods)
	writeTablesToFile(filename, [miscTable, periodTable, callsTables...])
end

function makeHospitalsStatsTables(periods::Vector{SimPeriodStats})::Vector{Table}
	numHospitals = length(periods[1].hospitals) # shorthand
	@assert(all(p -> length(p.hospitals) == numHospitals, periods))
	
	hospitalsTables = Table[]
	fnames = setdiff(fieldnames(HospitalStats), (:hospitalIndex,))
	getHospital(period::SimPeriodStats, hospitalIndex::Int) = hospitalIndex == 0 ? period.hospital : period.hospitals[hospitalIndex]
	for i = 0:numHospitals
		name = i == 0 ? "hospital" : "hospitals[$i]"
		hospitalTable = Table(name, vcat("periodIndex", collect(string.(fnames)));
			rows = [vcat(j, [getfield(getHospital(p,i), fname) for fname in fnames]) for (j,p) in enumerate(periods)])
		push!(hospitalsTables, hospitalTable)
	end
	
	return hospitalsTables
end

# for batches in a single simulation replication
function writeHospitalsStatsFile(filename::String, stats::SimStats)
	# shorthand
	periods = stats.periods
	numHospitals = length(periods[1].hospitals)
	
	miscTable = Table("miscData", ["mode", "numHospitals"]; rows = [["batch", numHospitals]])
	timestampsTable = simStatsTimestampsTable(stats)
	periodsTable = simStatsPeriodsTable(periods)
	hospitalsTables = makeHospitalsStatsTables(periods)
	writeTablesToFile(filename, [miscTable, timestampsTable, periodsTable, hospitalsTables...])
end

# for independent simulation replications
function writeHospitalsStatsFile(filename::String, periods::Vector{SimPeriodStats})
	numHospitals = length(periods[1].hospitals) # shorthand
	@assert(all(p -> length(p.hospitals) == numHospitals, periods))
	
	miscTable = Table("miscData", ["mode", "numHospitals"]; rows = [["replication", numHospitals]])
	periodTable = simStatsPeriodTable(periods)
	hospitalsTables = makeHospitalsStatsTables(periods)
	writeTablesToFile(filename, [miscTable, periodTable, hospitalsTables...])
end

function makeStationsStatsTables(periods::Vector{SimPeriodStats})::Vector{Table}
	numStations = length(periods[1].stations) # shorthand
	@assert(all(p -> length(p.stations) == numStations, periods))
	
	getStation(period::SimPeriodStats, stationIndex::Int) = stationIndex == 0 ? period.station : period.stations[stationIndex]
	stationsNumIdleAmbsTotalDurationTable = Table("stations_numIdleAmbsTotalDuration", vcat("periodIndex", "station", ["stations[$i]" for i = 1:numStations]);
		rows = [vcat(j, [string(getStation(p,i).numIdleAmbsTotalDuration) for i = 0:numStations]) for (j,p) in enumerate(periods)])
	return [stationsNumIdleAmbsTotalDurationTable]
end

# for batches in a single simulation replication
function writeStationsStatsFile(filename::String, stats::SimStats)
	# shorthand
	periods = stats.periods
	numStations = length(periods[1].stations)
	
	miscTable = Table("miscData", ["mode", "numStations"]; rows = [["batch", numStations]])
	timestampsTable = simStatsTimestampsTable(stats)
	periodsTable = simStatsPeriodsTable(periods)
	stationsTables = makeStationsStatsTables(periods)
	writeTablesToFile(filename, [miscTable, timestampsTable, periodsTable, stationsTables...])
end

# for independent simulation replications
function writeStationsStatsFile(filename::String, periods::Vector{SimPeriodStats})
	numStations = length(periods[1].stations) # shorthand
	@assert(all(p -> length(p.stations) == numStations, periods))
	
	miscTable = Table("miscData", ["mode", "numStations"]; rows = [["replication", numStations]])
	periodTable = simStatsPeriodTable(periods)
	stationsTables = makeStationsStatsTables(periods)
	writeTablesToFile(filename, [miscTable, periodTable, stationsTables...])
end

# stats dict values into table
function makeStatsDictTable(statsDict::Dict{String,Any})
	statsDictFlat = flatten(statsDict)
	@assert(all(v -> isa(v, MeanAndHalfWidth), values(statsDictFlat))) # check that statsDictFlat really is flat
	statsDictTable = Table("statsDict", ["stat", "mean", "halfWidth"];
		rows = [[k, statsDictFlat[k].mean, statsDictFlat[k].halfWidth] for k in sort(collect(keys(statsDictFlat)))])
	return statsDictTable
end

function writeStatsDictFile(filename::String, periods::Vector{SimPeriodStats}; conf::Float = 0.95)
	periodDuration = 0.0
	if !isempty(periods)
		periodDuration = periods[1].duration
		@assert(all(p -> isapprox(p.duration, periodDuration), periods))
	end
	statsDictFlat = statsDictFromPeriodStatsList(periods, conf = conf) |> flatten
	statsDictTable = makeStatsDictTable(statsDictFlat)
	miscTable = Table("misc", ["conf", "numPeriods", "periodDuration"], rows = [[conf, length(periods), periodDuration]])
	writeTablesToFile(filename, [miscTable, statsDictTable])
end
