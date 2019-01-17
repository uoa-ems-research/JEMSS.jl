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

# integer program formulation by Oddo Zhang, for move-up

# initialise data relevant to move up
function initZhangIp!(sim::Simulation;
	paramsFilename::String = "")
	
	# read parameters from file
	zid = sim.moveUpData.zhangIpData = readZhangIpParamsFile(paramsFilename)
	
	# check number of stations
	numStations = sim.numStations
	@assert(length(zid.marginalBenefits) == numStations)
	@assert(length(zid.stationCapacities) == numStations)
	
	# check station capacities
	for i = 1:numStations
		@assert(zid.stationCapacities[i] <= sim.stations[i].capacity)
	end
	
	for i = 1:numStations
		for j = 1:zid.stationCapacities[i]
			push!(zid.stationSlots, i)
			push!(zid.benefitSlots, zid.marginalBenefits[i][j])
		end
	end
	
	zid.marginalBenefitsDecreasing = all(i -> issorted(zid.marginalBenefits[i], lt=<=, rev=true), 1:numStations)
	if !zid.marginalBenefitsDecreasing
		p1 = vcat([findall(zid.stationSlots .== j)[1:end-1] for j = 1:numStations]...)
		p2 = vcat([findall(zid.stationSlots .== j)[2:end] for j = 1:numStations]...)
		zid.stationSlotsOrderPairs = hcat(p1, p2)
	end
	
	@assert(sim.travel.numSets == 1) # otherwise, need to be careful in calculating the regret travel time for move up of on-road ambulances
end

# determine move ups to make at current time
# returns list of ambulances to be moved, and list of their destinations (stations)
function zhangIpMoveUp(sim::Simulation)
	
	# shorthand names:
	zid = sim.moveUpData.zhangIpData
	stationSlots = zid.stationSlots
	benefitSlots = zid.benefitSlots
	ambulances = sim.ambulances
	stations = sim.stations
	numStations = sim.numStations
	
	# get currently movable ambulances, and at-hospital ambulances
	movableAmbs = filter(a -> isAmbAvailableForMoveUp(a), ambulances)
	atHospitalAmbs = filter(a -> a.status == ambAtHospital, ambulances)
	@assert(intersect(movableAmbs, atHospitalAmbs) == [])
	
	# let "move-up ambulances" be the ambulances that can be moved now, and those that can be moved later (currently at-hospital)
	moveUpAmbs = vcat(movableAmbs, atHospitalAmbs)
	numMoveUpAmbs = length(moveUpAmbs)
	
	# calculate travel time for each move-up ambulance to reach every station
	ambToStationTimes = Array{Float,2}(undef, numMoveUpAmbs, numStations)
	for i = 1:numMoveUpAmbs
		ambToStationTimes[i,:] = ambMoveUpTravelTimes!(sim, moveUpAmbs[i])
	end
	
	# calculate adjustedAmbToStationTimes from ambToStationTimes, according to Zhang's thesis
	# for at-station or newly-freed ambulances: adjusted time = original time
	# for on-road ambulances: adjusted time = original time * ( time spent on route + time to drive to station - time from route origin to station <= regret-travel-time threshold ? discount factor : 1 )
	# for at-hospital ambulances: adjusted time = original time + expected remaining transfer duration
	adjustedAmbToStationTimes = deepcopy(ambToStationTimes)
	ambIsNewlyIdle(a::Ambulance) = (a.status == ambGoingToStation && a.route.startTime == sim.time)
	travelMode = getTravelMode!(sim.travel, lowPriority, sim.time)
	for i = 1:numMoveUpAmbs
		ambulance = moveUpAmbs[i]
		if ambulance.status == ambIdleAtStation || ambIsNewlyIdle(ambulance)
			# no change, adjustedAmbToStationTimes[i,:] = ambToStationTimes[i,:]
		elseif ambulance.status == ambGoingToStation
			for j = 1:numStations
				# calculate "regret" travel time
				if ambulance.stationIndex == j
					regretTravelTime = 0.0
				else
					# get travel time from ambulance.route.startLoc to station j
					time1 = ambulance.route.startFNodeTime - ambulance.route.startTime
					pathTravelTime = shortestPathTravelTime(sim.net, travelMode.index, ambulance.route.startFNode, stations[j].nearestNodeIndex)
					time2 = offRoadTravelTime(travelMode, stations[j].nearestNodeDist)
					tob = time1 + pathTravelTime + time2 # on and off road travel times
					
					tox = sim.time - ambulance.route.startTime # time spent on current route
					txb = ambToStationTimes[i,j]
					regretTravelTime = tox + txb - tob
				end
				
				# apply discount if regret travel time is within threshold
				if regretTravelTime <= zid.regretTravelTimeThreshold
					adjustedAmbToStationTimes[i,j] *= zid.onRoadMoveUpDiscountFactor
				end
			end
		elseif ambulance.status == ambAtHospital
			adjustedAmbToStationTimes[i,:] += zid.expectedHospitalTransferDuration
		else
			error()
		end
	end
	
	travelCosts = adjustedAmbToStationTimes * zid.travelTimeCost
	
	######################
	# IP
	
	# shorthand variable names:
	a = numMoveUpAmbs
	ai = findin(moveUpAmbs, movableAmbs) # indices of movableAmbs in moveUpAmbs
	aj = findin(moveUpAmbs, atHospitalAmbs) # indices of atHospitalAmbs in moveUpAmbs
	s = numStations
	n = length(stationSlots) # = length(benefitSlots)
	
	# using JuMP
	model = Model()
	
	setsolver(model, GLPKSolverMIP(presolve=true))
	
	@variables(model, begin
		(x[i=1:a,j=1:s], Bin) # x[i,j] = 1 if ambulance moveUpAmbs[i] should be moved to station j
		(y[k=1:n]) # sum(y[stationSlots .== j]) = number of ambulances assigned to station j
	end)
	if zid.marginalBenefitsDecreasing
		# y should be naturally binary
		for k = 1:n setlowerbound(y[k], 0); setupperbound(y[k], 1) end
	else
		for k = 1:n setcategory(y[k], :Bin) end
	end
	
	@constraints(model, begin
		(movableAmbAtOneLocation[i=ai], sum(x[i,:]) == 1) # each movable ambulance must be assigned to one station
		(atHospitalAmbAtMostOneLocation[i=aj], sum(x[i,:]) <= 1) # each at-hospital ambulance may be assigned to zero or one station
		(stationAmbCounts[j=1:s], sum(x[:,j]) == sum(y[k] for k=findall(stationSlots .== j)))
		# (stationAmbCounts[j=1:s], sum(x[:,j]) >= sum(y[k] for k=findall(stationSlots .== j))) # should have same effect as "==" constraint (instead of ">="), but may be faster?
	end)
	
	if !zid.marginalBenefitsDecreasing
		# need to enforce station slot filling order
		p = zid.stationSlotsOrderPairs # shorthand
		@constraint(model, stationSlotsFillingOrder[k=1:size(p,1)], y[p[k,1]] >= y[p[k,2]])
	end
	
	@expressions(model, begin
		totalBenefitAtStations, sum(y .* benefitSlots)
		totalAmbTravelCosts, sum(x .* travelCosts)
	end)
	
	@objective(model, :Max, totalBenefitAtStations - totalAmbTravelCosts)
	
	# # testing: giving back fake results, for testing runtime without solving IP model
	# if true
		# return moveUpNull()
	# end
	
	solve(model)
	
	# extract solution
	sol = convert(Array{Bool,2}, round.(getvalue(x)))
	ambStations = [stations[findfirst(sol[i,:])] for i = ai] # only consider movableAmbs (ignore atHospitalAmbs)
	
	if checkMode
		@assert(all(sum(sol,2) .<= 1)) # each ambulance can be used in move up at most once
		
		# check that y values are ordered correctly
		stationSlotsFilled = convert(Vector{Bool}, round.(getvalue(y)))
		for j = 1:numStations
			@assert(issorted(stationSlotsFilled[stationSlots .== j], rev=true))
		end
	end
	
	return movableAmbs, ambStations
end
