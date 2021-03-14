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
# From thesis: "Simulation optimisation and Markov models for dynamic ambulance redeployment"

# initialise data relevant to move up
function initZhangIp!(sim::Simulation, zhangIpData::ZhangIpData)
	zid = sim.moveUpData.zhangIpData = zhangIpData # shorthand
	
	# check number of stations
	numStations = sim.numStations
	@assert(length(zid.marginalBenefits) == numStations)
	@assert(length(zid.stationCapacities) == numStations)
	
	# check station capacities
	for i = 1:numStations
		@assert(zid.stationCapacities[i] <= sim.stations[i].capacity)
	end
	
	zid.stationSlots = []
	zid.benefitSlots = []
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
function initZhangIp!(sim::Simulation; paramsFilename::String = "")
	zid = readZhangIpParamsFile(paramsFilename)
	initZhangIp!(sim, zid)
end

# determine move ups to make at current time
# returns list of ambulances to be moved, and list of their destinations (stations)
function zhangIpMoveUp(sim::Simulation)::Tuple{Vector{Ambulance}, Vector{Station}}
	
	# shorthand names:
	zid = sim.moveUpData.zhangIpData
	@unpack stationSlots, benefitSlots = zid
	@unpack ambulances, stations, numStations = sim
	
	# get currently movable ambulances, and at-hospital ambulances
	movableAmbs = filter(a -> isAmbMovable(a), ambulances)
	atHospitalAmbs = filter(a -> a.status == ambAtHospital, ambulances)
	@assert(!any(a -> isAmbMovable(a) && a.status == ambAtHospital, ambulances)) # @assert(isempty(intersect(movableAmbs, atHospitalAmbs))) is much slower
	if isempty(movableAmbs) return moveUpNull() end
	
	# let moveUpAmbs be the ambulances that can be moved now, and those that can be moved later (currently at-hospital)
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
	# for at-hospital ambulances: adjusted time = original time + expected remaining handover duration
	adjustedAmbToStationTimes = copy(ambToStationTimes)
	travelMode = getTravelMode!(sim.travel, lowPriority, sim.time)
	for i = 1:numMoveUpAmbs
		ambulance = moveUpAmbs[i]
		if ambulance.status == ambIdleAtStation || ambulance.status == ambFreeAfterCall
			# no change, adjustedAmbToStationTimes[i,:] = ambToStationTimes[i,:]
		elseif isGoingToStation(ambulance.status)
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
			adjustedAmbToStationTimes[i,:] .+= zid.expectedHospitalHandoverDuration
		else
			error()
		end
	end
	
	ambToStationCosts = adjustedAmbToStationTimes * zid.travelTimeCost
	
	# solve
	if zid.marginalBenefitsDecreasing
		# solve as assignment problem
		(movableAmbStations, atHospitalAmbStations) = solveZhangIpAssignmentProblem(stationSlots, benefitSlots, ambToStationCosts, length(movableAmbs), length(atHospitalAmbs))
	else
		# solve as an integer program
		(movableAmbStations, atHospitalAmbStations) = solveZhangIp(stationSlots, benefitSlots, ambToStationCosts, length(movableAmbs), length(atHospitalAmbs), zid.marginalBenefitsDecreasing;
			stationSlotsOrderPairs = zid.stationSlotsOrderPairs)
	end
	
	ambStations = stations[movableAmbStations]
	
	return movableAmbs, ambStations
end

# solve as an integer program
function solveZhangIp(stationSlots::Vector{Int}, benefitSlots::Vector{Float}, ambToStationCosts::Array{Float,2},
	numMovableAmbs::Int, numAtHospitalAmbs::Int,
	marginalBenefitsDecreasing::Bool; stationSlotsOrderPairs::Array{Int,2} = Array{Int,2}(undef,0,0))::Tuple{Vector{Int}, Vector{Int}}
	# ambToStationCosts[:,j] gives costs to move different ambulances to station j
	# ambToStationCosts[1:numMovableAmbs,:] are for ambs that can be moved now
	# ambToStationCosts[(1:numAtHospitalAmbs).+numMovableAmbs,:] are for ambs that are at hospital
	# stationSlotsOrderPairs is not needed if marginalBenefitsDecreasing = true
	
	(numMoveUpAmbs, numStations) = size(ambToStationCosts)
	if numMoveUpAmbs == 0 return (Int[], Int[]) end
	numStationSlots = length(stationSlots)
	@assert(numMoveUpAmbs == numMovableAmbs + numAtHospitalAmbs)
	@assert(numStationSlots >= numMoveUpAmbs)
	@assert(numStationSlots == length(benefitSlots))
	@assert(all(in(1:numStations), stationSlots))
	
	if checkMode
		if marginalBenefitsDecreasing
			for j = 1:numStations
				# check that benefit of ambulance k at station j is > benefit of ambulance k+1, otherwise marginalBenefitsDecreasing should be false
				@assert(issorted(benefitSlots[findall(stationSlots .== j)], lt=<=, rev=true))
			end
		end
	end
	
	# shorthand:
	a = numMoveUpAmbs
	ai = 1:numMovableAmbs # row indices of ambToStationCosts that are for movable ambs
	aj = (1:numAtHospitalAmbs).+numMovableAmbs # row indices of ambToStationCosts that are for at-hospital ambs
	s = numStations
	n = length(stationSlots) # = length(benefitSlots)
	
	# using JuMP
	model = Model()
	set_optimizer(model, with_optimizer(GLPK.Optimizer)) # faster without presolve; have not compared with other solvers
	# set_optimizer(model, with_optimizer(Cbc.Optimizer, logLevel=0))
	
	@variable(model, x[i=1:a,j=1:s], Bin)
	marginalBenefitsDecreasing ? @variable(model, 0 <= y[k=1:n] <= 1) : @variable(model, y[k=1:n], Bin)
	# x[i,j] = 1 if ambulance moveUpAmbs[i] should be moved to station j
	# sum(y[stationSlots .== j]) = number of ambulances assigned to station j
	
	@constraints(model, begin
		(movableAmbAtOneLocation[i=ai], sum(x[i,:]) == 1) # each movable ambulance must be assigned to one station
		(atHospitalAmbAtMostOneLocation[i=aj], sum(x[i,:]) <= 1) # each at-hospital ambulance may be assigned to zero or one station
		(stationAmbCounts[j=1:s], sum(x[:,j]) == sum(y[k] for k=findall(stationSlots .== j)))
		# (stationAmbCounts[j=1:s], sum(x[:,j]) >= sum(y[k] for k=findall(stationSlots .== j))) # should have same effect as "==" constraint (instead of ">="), but may be faster?
	end)
	
	if !marginalBenefitsDecreasing
		# need to enforce station slot filling order
		p = stationSlotsOrderPairs # shorthand
		@constraint(model, stationSlotsFillingOrder[k=1:size(p,1)], y[p[k,1]] >= y[p[k,2]])
	end
	
	@expressions(model, begin
		totalBenefitAtStations, sum(y .* benefitSlots)
		totalAmbTravelCosts, sum(x .* ambToStationCosts)
	end)
	
	# solve
	@objective(model, Max, totalBenefitAtStations - totalAmbTravelCosts)
	@stdout_silent optimize!(model)
	@assert(termination_status(model) == MOI.OPTIMAL)
	
	# get solution
	vals = Dict()
	vals[:x] = JuMP.value.(x)
	vals[:y] = JuMP.value.(y)
	
	# solution
	sol = convert(Array{Bool,2}, round.(vals[:x]))
	movableAmbStations = [findfirst(sol[i,:]) for i = ai]
	atHospitalAmbStations = [findfirst(sol[i,:]) for i = aj]
	atHospitalAmbStations = convert(Vector{Int}, replace(atHospitalAmbStations, nothing => nullIndex))
	
	if checkMode
		@assert(all(sum(sol, dims=2) .<= 1)) # each ambulance can be used in move up at most once
		
		# check that y values are ordered correctly
		stationSlotsFilled = convert(Vector{Bool}, round.(vals[:y]))
		for j = 1:numStations
			@assert(issorted(stationSlotsFilled[stationSlots .== j], rev=true))
		end
	end
	
	return movableAmbStations, atHospitalAmbStations # atHospitalAmbStations[i] == nullIndex if amb i (of those at hospital) is not moved
end

# solve as an assignment problem
# this assumes that the benefit of adding amb k to a station is > benefit of adding amb k+1
function solveZhangIpAssignmentProblem(stationSlots::Vector{Int}, benefitSlots::Vector{Float}, ambToStationCosts::Array{Float,2},
	numMovableAmbs::Int, numAtHospitalAmbs::Int)::Tuple{Vector{Int}, Vector{Int}}
	# ambToStationCosts[:,j] gives costs to move different ambulances to station j
	# ambToStationCosts[1:numMovableAmbs,:] are for ambs that can be moved now
	# ambToStationCosts[(1:numAtHospitalAmbs).+numMovableAmbs,:] are for ambs that are at hospital
	
	(numMoveUpAmbs, numStations) = size(ambToStationCosts)
	if numMoveUpAmbs == 0 return (Int[], Int[]) end
	@assert(numMoveUpAmbs == numMovableAmbs + numAtHospitalAmbs)
	numStationSlots = length(stationSlots)
	@assert(numStationSlots >= numMoveUpAmbs)
	@assert(numStationSlots == length(benefitSlots))
	@assert(all(in(1:numStations), stationSlots))
	
	if checkMode && false # skip this check, it is slow
		for j = 1:numStations
			# check that benefit of ambulance k at station j is > benefit of ambulance k+1, otherwise cannot solve as assignment problem
			@assert(issorted(benefitSlots[findall(stationSlots .== j)], lt=<=, rev=true))
		end
	end
	
	# shorthand
	a = numMoveUpAmbs
	ai = 1:numMovableAmbs # row indices of ambToStationCosts that are for movable ambs
	aj = (1:numAtHospitalAmbs).+numMovableAmbs # row indices of ambToStationCosts that are for at-hospital ambs
	n = numStationSlots + numAtHospitalAmbs # total number of slots, need one dummy per each amb at hospital
	si = 1:numStationSlots # column indices of ambToStationCosts that are for real stations
	sj = numStationSlots+1:n # column indices of ambToStationCosts that are for dummy stations
	
	# formulate assignment problem and solve
	@assert(a <= n) # needed for Hungarian algorithm
	weights = zeros(Float, a, n)
	for i = 1:a, j = si
		weights[i,j] = -(benefitSlots[j] - ambToStationCosts[i,stationSlots[j]]) # note that Hungarian algorithm looks for min cost assignment
	end
	weights[ai, sj] .= Inf # movable ambs cannot be assigned to dummy stations
	matching = Hungarian.munkres(weights) # returns a sparse matrix
	
	# extract solution
	I = findall(matching .== Hungarian.STAR)
	(ambIndices, stationSlotIndices) = (getindex.(I,1), getindex.(I,2))
	ambStations = zeros(Int,a) # will be the station (real or dummy) of each ambulance
	stationSlotsExtended = vcat(stationSlots, fill(nullIndex, numAtHospitalAmbs)) # add dummy station indices
	ambStations[ambIndices] = stationSlotsExtended[stationSlotIndices]
	movableAmbStations = ambStations[ai]
	atHospitalAmbStations = ambStations[aj]
	
	if checkMode
		@assert(!any(isequal(nullIndex), movableAmbStations)) # make sure that no movable ambs are assigned to dummy stations
		
		if false # skip this check, it is slow
			# check that station slots are filled in correct order
			stationSlotsFilled = fill(false, numStationSlots)
			for i in stationSlotIndices
				if i <= numStationSlots stationSlotsFilled[i] = true end
			end
			for j = 1:numStations
				@assert(issorted(stationSlotsFilled[findall(stationSlots .== j)], rev=true))
			end
		end
	end
	
	return movableAmbStations, atHospitalAmbStations # atHospitalAmbStations[i] == nullIndex if amb i (of those at hospital) is not moved
end
