# temp1 move-up

# initialise data relevant to move up
function initTemp1!(sim::Simulation;
	busyFraction::Float = 0.5, travelTimeCost::Float = 10.0, maxIdleAmbTravelTime::Float = 1.0, maxNumNearestStations::Int = 99)
	# shorthand names:
	tmp = sim.moveUpData.temp1Data
	numAmbs = sim.numAmbs
	stations = sim.stations
	numStations = sim.numStations
	
	# parameters:
	tmp.busyFraction = busyFraction
	tmp.travelTimeCost = travelTimeCost
	tmp.maxIdleAmbTravelTime = maxIdleAmbTravelTime # (days)
	tmp.maxNumNearestStations = maxNumNearestStations
	
	# create station pairs
	@assert(sim.travel.numSets == 1) # otherwise, would need to create station pairs for each travel set
	travelMode = getTravelMode!(sim.travel, lowPriority, sim.time)
	stationPairs = createStationPairs(sim, travelMode;
		maxPairsPerStation = 5, maxPairSeparation = 50/(24*60))
	
	# for easy indexing:
	stationSingles = Vector{Vector{Int}}(0)
	for i = 1:numStations
		push!(stationSingles, [i,i])
	end
	
	# testing / temporary:
	benefit = Array{Vector{Float},2}(numStations,numStations)
	for i = 1:numStations
		n = stations[i].capacity
		benefit[i,i] = 1 - busyFraction.^[1:n;]
	end
	for (i,j) in stationPairs
		n = stations[i].capacity + stations[j].capacity
		benefit[i,j] = 0.2*(1 - busyFraction.^[1:n;])
	end
	marginalBenefit = deepcopy(benefit)
	for (i,j) in [stationSingles; stationPairs]
		marginalBenefit[i,j][1] = benefit[i,j][1]
		for k = 2:length(marginalBenefit[i,j])
			marginalBenefit[i,j][k] = benefit[i,j][k] - benefit[i,j][k-1]
		end
	end
	
	tmp.benefit = benefit
	tmp.marginalBenefit = marginalBenefit
	tmp.stationPairs = stationPairs
end

# determine move ups to make at current time
# returns list of ambulances to be moved, and list of their destinations (stations)
function temp1MoveUp(sim::Simulation)

	# shorthand names:
	tmp = sim.moveUpData.temp1Data
	stationPairs = tmp.stationPairs
	travelTimeCost = tmp.travelTimeCost
	maxIdleAmbTravelTime = tmp.maxIdleAmbTravelTime
	maxNumNearestStations = tmp.maxNumNearestStations
	marginalBenefit = tmp.marginalBenefit
	ambulances = sim.ambulances
	stations = sim.stations
	
	numAmbs = sim.numAmbs
	numStations = sim.numStations

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
		for k = 1:numSlots
			push!(stationSlots, j)
			push!(benefitSlots, marginalBenefit[j,j][k])
		end
	end
	
	# counting number of ambulances at station pairs
	stationPairSlots = Vector{Int}(0)
	benefitPairSlots = Vector{Float}(0)
	mapmbz = Array{Vector{Int},2}(size(marginalBenefit)) # for mapping marginalBenefit to stationPairSlots
	for (p, (i,j)) in enumerate(stationPairs)
		numPairSlots = sum(stationSlots .== i) + sum(stationSlots .== j)
		mapmbz[i,j] = nullIndex * ones(Int, length(marginalBenefit[i,j]))
		for k = 1:numPairSlots
			push!(stationPairSlots, p)
			push!(benefitPairSlots, marginalBenefit[i,j][k])
			mapmbz[i,j][k] = length(stationPairSlots)
		end
	end

	######################
	# IP

	# shorthand variable names:
	a = numMovableAmbs
	s = numStations
	m = length(stationList)
	n = length(stationSlots) # <= m
	p = length(stationPairs)
	q = length(stationPairSlots)

	# using JuMP
	model = Model()
	
	setsolver(model, GLPKSolverMIP(presolve=true))

	@variables(model, begin
		(x[k=1:m], Bin) # x[k] = 1 if ambulance: ambList[k] should be moved to station: stationList[k]
		(0 <= y[k=1:n] <= 1) # y should be naturally binary; sum(y[stationSlots .== k]) = number of ambulances assigned to stations[k]
		(0 <= z[k=1:q] <= 1) # z should be naturally binary; sum(z[stationPairSlots .== k]) = number of ambulances assigned to station pair k
	end)

	@constraints(model, begin
		(ambAtOneLocation[i=1:a], sum(x[k] for k=find(ambList .== i)) == 1) # each ambulance must be assigned to one station
		(stationAmbCounts[j=1:s], sum(x[k] for k=find(stationList .== j)) == sum(y[k] for k=find(stationSlots .== j)))
		(stationPairAmbCounts[l=1:p], sum(x[k] for k=find(stationList .== stationPairs[l][1])) + sum(x[k] for k=find(stationList .== stationPairs[l][2])) == sum(z[k] for k=find(stationPairSlots .== l)))
	end)

	@expressions(model, begin
		totalAmbTravelCosts, sum(x[k] * travelCostList[k] for k=1:m)
		totalBenefitAtStations, sum(y[k] * benefitSlots[k] for k=1:n)
		totalBenefitAtStationPairs, sum(z[k] * benefitPairSlots[k] for k=1:q)
	end)

	@objective(model, :Max, totalBenefitAtStations + totalBenefitAtStationPairs - totalAmbTravelCosts)

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
			@assert(issorted(stationSlotsFilled[stationSlots .== i], rev=true))
		end
		
		# check that z values are ordered correctly
		stationPairSlotsFilled = convert(Vector{Bool}, round.(getvalue(z)))
		for (i,j) in stationPairs
			f(k) = mapmbz[i,j][k] == nullIndex ? 0 : stationPairSlotsFilled[mapmbz[i,j][k]]
			n = length(marginalBenefit[i,j])
			for k = 2:n
				@assert(f(k) <= f(k-1))
			end
		end
	end

	return movableAmbs, ambStations
end
