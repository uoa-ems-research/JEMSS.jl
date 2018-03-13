# integer program by Oddo Zhang, for move-up

# todo: should allow move up of soon to be idle ambulances, as done by Zhang

# initialise data relevant to move up
function initZhangIp!(sim::Simulation;
	busyFraction::Float = 0.5, travelTimeCost::Float = 10.0, maxIdleAmbTravelTime::Float = 1.0, maxNumNearestStations::Int = 99)
	# shorthand names:
	zid = sim.moveUpData.zhangIpData
	numAmbs = length(sim.ambulances)
	
	# parameters:
	zid.busyFraction = busyFraction
	zid.travelTimeCost = travelTimeCost
	zid.maxIdleAmbTravelTime = maxIdleAmbTravelTime # (days)
	zid.maxNumNearestStations = maxNumNearestStations
	
	zid.marginalBenefit  = (busyFraction.^[0:numAmbs-1;])*(1-busyFraction)
end

# determine move ups to make at current time
# returns list of ambulances to be moved, and list of their destinations (stations)
function zhangIpMoveUp(sim::Simulation)

	# shorthand names:
	zid = sim.moveUpData.zhangIpData
	travelTimeCost = zid.travelTimeCost
	maxIdleAmbTravelTime = zid.maxIdleAmbTravelTime
	maxNumNearestStations = zid.maxNumNearestStations
	marginalBenefit = zid.marginalBenefit
	ambulances = sim.ambulances
	stations = sim.stations
	
	numAmbs = length(ambulances)
	numStations = length(stations)

	# get movable ambulances (movableAmbs)
	ambMovable = Vector{Bool}(numAmbs) # ambMovable[i] = true if ambulances[i] can move-up
	for i = 1:numAmbs
		ambMovable[i] = isAmbAvailableForMoveUp(ambulances[i])
	end
	movableAmbs = ambulances[ambMovable]
	numMovableAmbs = length(movableAmbs)

	# calculate travel time for each available ambulance to reach every station
	ambToStationTimes = Array{Float,2}(numMovableAmbs, numStations)
	for i = 1:numMovableAmbs
		ambToStationTimes[i,:] = ambMoveUpTravelTimes!(sim, movableAmbs[i])
	end
	# edit travel times so that an ambulance travelling to its own station has no travel time
	for i = 1:numMovableAmbs
		j = movableAmbs[i].stationIndex
		ambToStationTimes[i,j] = 0.0
	end

	# restrict which stations each ambulance can be moved to
	# ambMovableToStation[i,j] = true if movableAmbs[i] can be moved to stations[j]; false otherwise
	ambMovableToStation = Array{Bool,2}(numMovableAmbs, numStations)
	ambMovableToStation[:,:] = true

	# limit ambulance move-up to nearest stations
	numNearestStations = min(maxNumNearestStations, numStations)
	sortedTimes = sort(ambToStationTimes, 2) # sortedTimes[i,:] gives travel times for ith movable ambulance to all other stations, sorted
	ambMovableToStation = (ambMovableToStation .* (ambToStationTimes .<= sortedTimes[:, numNearestStations]))

	# for ambulances idle at station, limit travel time
	for i = 1:numMovableAmbs
		ambulance = movableAmbs[i]
		if ambulance.status == ambIdleAtStation
			ambMovableToStation[i,:] = (ambMovableToStation[i,:] .* (ambToStationTimes[i,:] .<= maxIdleAmbTravelTime))
		end
	end

	# useful lists for IP
	(ambList, stationList) = findn(ambMovableToStation)
	# ambList and stationList together have all the information of ambMovableToStation:
	# - movableAmbs[i] can move to stations stationList[ambList .== i]
	# - stations[j] can have any of the ambulances in ambList[stationList .== j]
	m = length(ambList) # number of variables needed for assignment of ambulances to stations
	travelCostList = Vector{Float}(m)
	for k = 1:m
		travelCostList[k] = ambToStationTimes[ambList[k], stationList[k]] * travelTimeCost
	end

	# counting number of ambulances at each station
	stationSlots = Vector{Int}(0)
	benefitSlots = Vector{Float}(0)
	for j = 1:numStations
		numSlots = min(stations[j].capacity, sum(stationList .== j))
		for i = 1:numSlots
			push!(stationSlots, j)
			push!(benefitSlots, marginalBenefit[i])
		end
	end

	######################
	# IP

	# shorthand variable names:
	a = numMovableAmbs
	s = numStations
	m = length(stationList)
	n = length(stationSlots) # <= m

	# using JuMP
	model = Model()
	
	setsolver(model, GLPKSolverMIP(presolve=true))
	
	@variables(model, begin
		(x[k=1:m], Bin) # x[k] = 1 if ambulance: ambList[k] should be moved to station: stationList[k]
		(0 <= y[l=1:n] <= 1) # y should be naturally binary; sum(y[stationSlots .== k]) = number of ambulances assigned to stations[k]
	end)

	@constraints(model, begin
		(ambAtOneLocation[i=1:a], sum(x[k] for k=find(ambList .== i)) == 1) # each ambulance must be assigned to one station
		(stationAmbCounts[j=1:s], sum(x[k] for k=find(stationList .== j)) == sum(y[l] for l=find(stationSlots .== j)))
	end)

	@expressions(model, begin
		totalAmbTravelCosts, sum(x[k] * travelCostList[k] for k=1:m)
		totalBenefitAtStations, sum(y[l] * benefitSlots[l] for l=1:n)
	end)

	@objective(model, :Max, totalBenefitAtStations - totalAmbTravelCosts)

	# # testing: giving back fake results, for testing runtime without solving IP model
	# if true
		# return moveUpNull()
	# end

	solve(model)

	# extract solution
	sol = convert(Vector{Bool}, round.(getvalue(x)))
	ambStations = Vector{Station}(numMovableAmbs)
	for k = 1:m
		if sol[k] == 1
			ambStations[ambList[k]] = stations[stationList[k]]
		end
	end
	
	if checkMode
		# check that y values are ordered correctly
		stationSlotsFilled = convert(Vector{Bool}, round.(getvalue(y)))
		for i = 1:numStations
			assert(issorted(stationSlotsFilled[stationSlots .== i], rev=true))
		end
	end

	return movableAmbs, ambStations
end
