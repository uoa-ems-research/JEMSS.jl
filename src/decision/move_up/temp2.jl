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

# temp2 move-up

# initialise data relevant to move up
function initTemp2!(sim::Simulation;
	busyFraction::Float = 0.5, travelTimeCost::Float = 10.0, maxIdleAmbTravelTime::Float = 1.0, maxNumNearestStations::Int = 99)
	
	pkgVersions["JuMP"] >= v"0.19" && @warn("Move up module 'temp2' will not work with JuMP v0.19 and above.")
	
	# shorthand names:
	tmp = sim.moveUpData.temp2Data
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
		maxPairsPerStation = 10, maxPairSeparation = 60/(24*60))
	
	# for easy indexing:
	stationSingles = Vector{Vector{Int}}()
	for i = 1:numStations
		push!(stationSingles, [i,i])
	end
	
	# testing / temporary:
	marginalBenefit = Array{Array{Float,2},2}(undef,numStations,numStations)
	for (i,j) in [stationSingles; stationPairs]
		m = stations[i].capacity
		n = stations[j].capacity
		mb = Array{Float,2}(undef,m,n)
		mb[1:m,1] = (1-busyFraction)*busyFraction.^[0:m-1;]
		mb[1,1:n] = (1-busyFraction)*busyFraction.^[0:n-1;]
		if i != j
			mb *= 0.2
		end
		for k1 = 2:m, k2 = 2:n
			mb[k1,k2] = min(mb[k1,k2-1], mb[k1-1,k2]) * 0.5
		end
		marginalBenefit[i,j] = mb
	end
	benefit = deepcopy(marginalBenefit)
	for (i,j) in [stationSingles; stationPairs]
		m = stations[i].capacity
		n = stations[j].capacity
		for k1 = 1:m, k2 = 1:n
			benefit[i,j][k1,k2] = sum(marginalBenefit[i,j][1:k1,1:k2])
		end
	end
	
	# note that for 'benefit' and 'marginalBenefit', if the station indices are the same, only the first row/column should be used
	tmp.benefit = benefit
	tmp.marginalBenefit = marginalBenefit
	tmp.stationPairs = stationPairs
end

# determine move ups to make at current time
# returns list of ambulances to be moved, and list of their destinations (stations)
function temp2MoveUp(sim::Simulation)
	
	# shorthand names:
	tmp = sim.moveUpData.temp2Data
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
	ambMovable = Vector{Bool}(undef, numAmbs) # ambMovable[i] = true if ambulances[i] can move-up
	for i = 1:numAmbs
		ambMovable[i] = isAmbAvailableForMoveUp(ambulances[i])
	end
	movableAmbs = ambulances[ambMovable]
	numMovableAmbs = length(movableAmbs)
	
	# calculate travel time for each available ambulance to reach every station
	ambToStationTimes = Array{Float,2}(undef, numMovableAmbs, numStations)
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
	ambMovableToStation = Array{Bool,2}(undef, numMovableAmbs, numStations)
	ambMovableToStation[:,:] .= true
	
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
	I = findall(ambMovableToStation)
	(ambList, stationList) = (getindex.(I, 1), getindex.(I, 2))
	# ambList and stationList together have all the information of ambMovableToStation:
	# - movableAmbs[i] can move to stations stationList[ambList .== i]
	# - stations[j] can have any of the ambulances in ambList[stationList .== j]
	m = length(ambList) # number of variables needed for assignment of ambulances to stations
	travelCostList = Vector{Float}(undef, m)
	for k = 1:m
		travelCostList[k] = ambToStationTimes[ambList[k], stationList[k]] * travelTimeCost
	end
	
	# counting number of ambulances at each station
	stationSlots = Vector{Int}()
	benefitSlots = Vector{Float}()
	for j = 1:numStations
		numSlots = min(stations[j].capacity, sum(stationList .== j))
		for k = 1:numSlots
			push!(stationSlots, j)
			push!(benefitSlots, marginalBenefit[j,j][1,k])
		end
	end
	
	# counting product of number of ambulances at station pairs
	stationPairSlots = Vector{Int}()
	benefitPairSlots = Vector{Float}()
	mapzy = Vector{Vector{Int}}() # for mapping stationPairSlots to pair of stationSlots
	mapmbz = Array{Array{Int,2},2}(undef, size(marginalBenefit)) # for mapping marginalBenefit to stationPairSlots
	for (p, (i,j)) in enumerate(stationPairs)
		ssi = findall(stationSlots .== i)
		ssj = findall(stationSlots .== j)
		mapmbz[i,j] = nullIndex * ones(Int, size(marginalBenefit[i,j]))
		for k1 = 1:length(ssi), k2 = 1:length(ssj)
			push!(stationPairSlots, p)
			push!(benefitPairSlots, marginalBenefit[i,j][k1,k2])
			push!(mapzy, [ssi[k1], ssj[k2]])
			mapmbz[i,j][k1,k2] = length(stationPairSlots)
		end
	end
	
	# testing: calculate lower bound on total benefit
	numAmbsAtStation = zeros(Int,numStations)
	for i = 1:numMovableAmbs
		numAmbsAtStation[movableAmbs[i].stationIndex] += 1
	end
	totalBenefitLowerBound = 0.0
	for i = 1:numStations
		totalBenefitLowerBound += sum(marginalBenefit[i,i][1,1:numAmbsAtStation[i]])
	end
	for (i,j) in stationPairs
		totalBenefitLowerBound += sum(sum(marginalBenefit[i,j][1:numAmbsAtStation[i],1:numAmbsAtStation[j]]))
	end
	
	######################
	# IP
	
	# shorthand variable names:
	a = numMovableAmbs
	s = numStations
	m = length(stationList)
	n = length(stationSlots) # <= m
	q = length(stationPairSlots) # = length(benefitPairSlots)
	
	# using JuMP
	model = Model()
	
	setsolver(model, GLPKSolverMIP(presolve=true))
	
	@variables(model, begin
		(x[k=1:m], Bin) # x[k] = 1 if ambulance: ambList[k] should be moved to station: stationList[k]
		(0 <= y[k=1:n] <= 1) # y should be naturally binary; sum(y[stationSlots .== k]) = number of ambulances assigned to stations[k]
		(0 <= z[k=1:q] <= 1) # z should be naturally binary
	end)
	
	@expressions(model, begin
		totalAmbTravelCosts, sum(x[k] * travelCostList[k] for k=1:m)
		totalBenefitAtStations, sum(y[k] * benefitSlots[k] for k=1:n)
		totalBenefitAtStationPairs, sum(z[k] * benefitPairSlots[k] for k=1:q)
		# totalBenefit, totalBenefitAtStations + totalBenefitAtStationPairs
	end)
	
	@constraints(model, begin
		(ambAtOneLocation[i=1:a], sum(x[k] for k=findall(ambList .== i)) == 1) # each ambulance must be assigned to one station
		(stationAmbCounts[j=1:s], sum(x[k] for k=findall(stationList .== j)) == sum(y[k] for k=findall(stationSlots .== j)))
		(stationPairAmbCounts[k=1:q, l=1:2], z[k] <= y[mapzy[k][l]])
		# (totalBenefit >= totalBenefitLowerBound)
	end)
	
	@objective(model, :Max, totalBenefitAtStations + totalBenefitAtStationPairs - totalAmbTravelCosts)
	
	# # testing: giving back fake results, for testing runtime without solving IP model
	# if true
		# return moveUpNull()
	# end
	
	solve(model)
	
	# extract solution
	sol = convert(Vector{Bool}, round.(getvalue(x)))
	ambStations = Vector{Station}(undef, numMovableAmbs)
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
			f(k1,k2) = mapmbz[i,j][k1,k2] == nullIndex ? 0 : stationPairSlotsFilled[mapmbz[i,j][k1,k2]]
			(m,n) = size(marginalBenefit[i,j])
			for k1 = 1:m, k2 = 1:n
				if k1 > 1
					@assert(f(k1,k2) <= f(k1-1,k2))
				end
				if k2 > 1
					@assert(f(k1,k2) <= f(k1,k2-1))
				end
			end
		end
	end
	
	return movableAmbs, ambStations
end
