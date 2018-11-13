# compliance table for move up

# initialise data relevant to move up
function initCompTable!(sim::Simulation, compTableFilename::String)
	# shorthand names:
	numAmbs = sim.numAmbs
	numStations = sim.numStations
	ctd = sim.moveUpData.compTableData
	
	# read compliance table
	ctd.compTable = readCompTableFile(compTableFilename)
	@assert(size(ctd.compTable) == (numAmbs, numStations))
	
	# compliance table should not violate station capacities
	for j = 1:numStations
		@assert(all(ctd.compTable[:,j] .<= sim.stations[j].capacity))
	end
	
	# set array sizes
	ctd.ambMovable = Vector{Bool}(numAmbs)
end

function compTableMoveUp(sim::Simulation)
	@assert(sim.moveUpData.useMoveUp)
	
	# shorthand:
	ambulances = sim.ambulances
	stations = sim.stations
	numAmbs = sim.numAmbs
	numStations = sim.numStations
	ctd = sim.moveUpData.compTableData
	compTable = ctd.compTable
	ambMovable = ctd.ambMovable

	# get movable ambulances (movableAmbs)
	for i = 1:numAmbs
		ambMovable[i] = isAmbAvailableForMoveUp(ambulances[i])
	end
	movableAmbs = ambulances[ambMovable]
	numMovableAmbs = length(movableAmbs)
	
	if numMovableAmbs == 0
		return [], []
	end
	
	# calculate travel time for each available ambulance to reach every station
	ambToStationTimes = Array{Float,2}(numMovableAmbs, numStations) # ambToStationTimes[i,j] = time for ambulance i to travel to station j
	for i = 1:numMovableAmbs
		ambToStationTimes[i,:] = ambMoveUpTravelTimes!(sim, movableAmbs[i])
	end
	# edit travel times so that an ambulance travelling to its own station has no travel time
	for i = 1:numMovableAmbs
		j = movableAmbs[i].stationIndex
		ambToStationTimes[i,j] = 0.0
	end

	######################
	# IP

	# shorthand variable names:
	a = numMovableAmbs
	s = numStations
	
	# solve IP to get available ambulances to match compliance table, while minimising driving time
	model = Model()
	
	setsolver(model, GLPKSolverMIP(presolve=true))
	
	@variables(model, begin
		(x[i=1:a,j=1:s], Bin) # x[i,j] = 1 if ambulance: movableAmbs[i] should be moved to station: stations[j]
		# maxTravelTime # if minimising max travel time of all ambs, instead of total
	end)
	
	@constraints(model, begin
		(followCompTable[j=1:s], sum(x[i,j] for i=1:a) == compTable[numMovableAmbs, j])
		(ambAtOneLocation[i=1:a], sum(x[i,j] for j=1:s) == 1) # each ambulance must be assigned to one station
		# maxTravelTimeBound[i=1:a], sum(x[i,j] * ambToStationTimes[i,j] for j=1:s) <= maxTravelTime # max travel time of all ambulances
	end)
	
	@expressions(model, begin
		totalAmbTravelTime, sum(x[i,j] * ambToStationTimes[i,j] for i=1:a, j=1:s)
	end)
	
	@objective(model, :Min, totalAmbTravelTime)
	# @objective(model, :Min, maxTravelTime)
	
	# # testing: giving back fake results, for testing runtime without solving IP model
	# if true
		# return moveUpNull()
	# end
	
	solve(model)
	
	# extract solution
	sol = convert(Array{Bool,2}, round.(getvalue(x)))
	ambStations = Vector{Station}(numMovableAmbs)
	for i = 1:a, j = 1:s
		if sol[i,j] == 1
			ambStations[i] = stations[j]
		end
	end
	
	if checkMode
		# check that followCompTable constraint was met
		for j = 1:numStations
			@assert(sum(sol[:,j]) == compTable[numMovableAmbs,j])
		end
		# check that ambAtOneLocation constraint was met
		for i = 1:numMovableAmbs
			@assert(sum(sol[i,:]) == 1)
		end
	end
	
	return movableAmbs, ambStations
end
