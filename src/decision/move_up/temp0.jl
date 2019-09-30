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

# temp0 move up

# initialise data relevant to move up
function initTemp0!(sim::Simulation;
	busyFraction::Float = 0.5, travelTimeCost::Float = 10.0, maxIdleAmbTravelTime::Float = 1.0, maxNumNearestStations::Int = 99)
	
	pkgVersions["JuMP"] >= v"0.19" && @warn("Move up module 'temp0' will not work with JuMP v0.19 and above.")
	
	# shorthand names:
	tmp = sim.moveUpData.temp0Data
	numAmbs = sim.numAmbs
	
	# parameters:
	tmp.busyFraction = busyFraction
	tmp.travelTimeCost = travelTimeCost
	tmp.maxIdleAmbTravelTime = maxIdleAmbTravelTime # (days)
	tmp.maxNumNearestStations = maxNumNearestStations
	
	tmp.marginalBenefit = (busyFraction.^[0:numAmbs-1;])*(1-busyFraction)
end

# determine move ups to make at current time
# returns list of ambulances to be moved, and list of their destinations (stations)
function temp0MoveUp(sim::Simulation)
	
	# shorthand names:
	tmp = sim.moveUpData.temp0Data
	@unpack travelTimeCost, maxIdleAmbTravelTime, maxNumNearestStations, marginalBenefit = tmp
	@unpack ambulances, stations, numAmbs, numStations = sim
	
	# get movable ambulances (movableAmbs)
	ambMovable = Vector{Bool}(undef, numAmbs) # ambMovable[i] = true if ambulances[i] can move-up
	for i = 1:numAmbs
		ambMovable[i] = isAmbMovable(ambulances[i])
	end
	movableAmbs = ambulances[ambMovable]
	numMovableAmbs = length(movableAmbs)
	
	# calculate travel time for each movable ambulance to reach every station
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
		(ambAtOneLocation[i=1:a], sum(x[k] for k=findall(ambList .== i)) == 1) # each ambulance must be assigned to one station
		(stationAmbCounts[j=1:s], sum(x[k] for k=findall(stationList .== j)) == sum(y[l] for l=findall(stationSlots .== j)))
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
	end
	
	return movableAmbs, ambStations
end
