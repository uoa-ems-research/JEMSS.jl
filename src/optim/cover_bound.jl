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

# Cover bound - an upper bound on the fraction of responses in time for optimal move-up.
# From paper: "A bound on the performance of an optimal ambulance redeployment policy".

# mutates: coverBound, coverBoundSim
function calcCoverBound!(sim::Simulation;
		coverBound::Union{CoverBound,Nothing} = nothing,
		coverBoundSim::Union{CoverBoundSim,Nothing} = nothing, # need this if coverBound.sim is not set
		ambBusyDurationsToSample::Vector{Float} = Float[],
		queuedDurationsToSample::Vector{Float} = Float[],
		doPrint::Bool = false)
	
	@assert(ambBusyDurationsToSample != [] && all(ambBusyDurationsToSample .> 0))
	@assert(issorted(ambBusyDurationsToSample))
	
	if coverBound == nothing
		coverBound = initCoverBound(sim, doPrint = doPrint)
		@assert(coverBoundSim != nothing)
		coverBound.sim = coverBoundSim
	end
	cb = coverBound # shorthand
	cbs = coverBoundSim = coverBound.sim # shorthand
	
	if cb.ambBusyDurationsToSample != ambBusyDurationsToSample || length(cb.ambBusyDurationsToSample) != size(cb.ambBusyDurationProbUpperBounds, 1)
		cb.ambBusyDurationsToSample = ambBusyDurationsToSample
		cb.ambBusyDurationProbUpperBounds = calcAmbBusyDurationProbUpperBounds(sim, cb, doPrint = doPrint)
	else
		# user should ensure that ambBusyDurationProbUpperBounds are for ambBusyDurationsToSample
		doPrint && println("Using existing coverBound.ambBusyDurationProbUpperBounds.")
	end
	
	calcAmbBusyDurationLowerBoundDistrs!(cb) # create amb busy duration distributions
	
	cb.accountForQueuedDurations = !isempty(queuedDurationsToSample)
	if cb.accountForQueuedDurations
		if cb.queuedDurationsToSample != queuedDurationsToSample || length(cb.queuedDurationsToSample) != length(cb.queuedDurationsMaxCoverageFrac)
			cb.queuedDurationsToSample = queuedDurationsToSample
			cb.queuedDurationsMaxCoverageFrac = calcQueuedDurationsMaxCoverageFrac(sim, cb, doPrint = doPrint)
		else
			# user should ensure that queuedDurationsMaxCoverageFrac are for queuedDurationsToSample
			doPrint && println("Using existing coverBound.queuedDurationsMaxCoverageFrac.")
		end
	end
	
	simulateCoverBound!(cb)
	
	return coverBound
end

# set kwarg pMedianRelax = true if using linear relaxation of p-median problem, false for binary
function initCoverBound(sim::Simulation; pMedianRelax::Bool = true, doPrint::Bool = false)
	
	coverBound = CoverBound()
	
	sim.demand.initialised || initDemand!(sim)
	sim.demandCoverage.initialised || initDemandCoverage!(sim)
	
	@assert(sim.demand.numSets == 1) # otherwise would have to calculate cover bound for each demand set
	@assert(sim.travel.numSets == 1) # otherwise would have to calculate cover bound for each travel set
	
	genConfig = readGenConfig(sim.inputFiles["callGenConfig"].path) # used for service time related distributions
	# will assume that all demand priorities have the same distributions
	
	coverBound.modeLookup = fill(nullIndex, sim.travel.numModes)
	
	numModes = 0
	for demandPriority in priorities
		responseTravelPriority = sim.responseTravelPriorities[demandPriority]
		travelMode = getTravelMode!(sim.travel, responseTravelPriority, sim.startTime) # assumes sim.travel.numSets == 1
		
		# check if there is already another coverBoundMode with the same response travel priority
		if coverBound.modeLookup[travelMode.index] != nullIndex
			continue # can skip this demandPriority, it will give the same coverBoundMode as another
		end
		
		coverBoundMode = initCoverBoundMode(sim, travelMode)
		numModes += 1
		coverBoundMode.index = numModes
		coverBound.modeLookup[travelMode.index] = numModes
		push!(coverBound.modes, coverBoundMode)
		
		setCoverBoundModeDistrs!(coverBoundMode, genConfig)
		calcApproxDistrs!(coverBoundMode) # pMax::Float = 0.999, dt::Float = 1/60/60/24
	end
	
	# will use deterministic dispatch delay, to simplify implementation
	dispatchDelayDistr = genConfig.dispatchDelayDistrRng.d
	if var(dispatchDelayDistr) != 0
		@warn("Have changed input dispatch delay distribution, as this implementation of the cover bound assumes deterministic dispatch delay.")
		dispatchDelayDistr = Normal(mean(dispatchDelayDistr), 0)
	end
	coverBound.dispatchDelay = mean(dispatchDelayDistr)
	@assert(coverBound.dispatchDelay >= 0)
	
	coverBound.numAmbsMaxCoverageFrac = calcNumAmbsMaxCoverageFrac(sim)
	
	return coverBound
end

function initCoverBoundMode(sim::Simulation, responseTravelMode::TravelMode)
	
	# shorthand
	currentTime = sim.startTime
	hospitalTravelPriority = lowPriority # default
	
	# shorthand:
	@unpack stations, numStations = sim
	numNodes = length(sim.net.fGraph.nodes)
	points = sim.demandCoverage.points
	numPoints = length(points)
	nodesPoints = sim.demandCoverage.nodesPoints
	
	# calculate travel time from every station to every node and point
	travelMode = responseTravelMode # shorthand
	stationsToNodesTimes = zeros(Float, numStations, numNodes)
	stationsToPointsTimes = zeros(Float, numStations, numPoints)
	for i = 1:numStations
		station = stations[i] # shorthand
		
		# get travel time from station to nearest node
		(node1, dist1) = (station.nearestNodeIndex, station.nearestNodeDist)
		time1 = offRoadTravelTime(travelMode, dist1) # time to reach nearest node from station
		
		for j = 1:numNodes
			pathTravelTime = shortestPathTravelTime(sim.net, travelMode.index, node1, j)
			stationsToNodesTimes[i,j] = time1 + pathTravelTime
			for k in nodesPoints[j] # indices of points for which node j is the nearest node
				point = points[k] # shorthand
				@assert(point.nearestNodeIndex == j)
				time2 = offRoadTravelTime(travelMode, point.nearestNodeDist) # time to reach point from nearest node
				stationsToPointsTimes[i,k] = time1 + pathTravelTime + time2
			end
		end
	end
	
	# group nodes based on the order in which each station can reach them (closest to furthest)
	nodeSetsDict = Dict{Vector{Int}, Vector{Int}}() # nodeSetsDict[perm] gives the node indices that can be reached by stations in order stations[perm]
	for j = 1:numNodes
		stationsToNodeTimes = view(stationsToNodesTimes, :, j)
		perm = sortperm(stationsToNodeTimes)
		v = get!(nodeSetsDict, perm, Int[])
		push!(v, j)
	end
	nodeSets = collect(values(nodeSetsDict))
	nodeSetsStationList = collect(keys(nodeSetsDict)) # nodeSetsStationList[i] = list of station indices in order (closest to furthest) that they cover nodeSets[i]
	
	# calculate travel time from every demand point to the closest hospital
	travelMode = getTravelMode!(sim.travel, hospitalTravelPriority, currentTime) # assumes travel.numSets == 1
	pointsToHospitalTimes = zeros(Float, numPoints)
	hospitalToNearestNodeTime = [offRoadTravelTime(travelMode, h.nearestNodeDist) for h in sim.hospitals]
	for j = 1:numNodes
		hospitalIndex = travelMode.fNetTravel.fNodeNearestHospitalIndex[j]
		node2 = sim.hospitals[hospitalIndex].nearestNodeIndex
		pathTravelTime = shortestPathTravelTime(sim.net, travelMode.index, j, node2)
		time2 = hospitalToNearestNodeTime[hospitalIndex]
		for k in nodesPoints[j] # indices of points for which node j is the nearest node
			point = points[k] # shorthand
			time1 = offRoadTravelTime(travelMode, point.nearestNodeDist) # time to reach point from nearest node
			@assert(pointsToHospitalTimes[k] == 0) # value should not yet be initialised
			pointsToHospitalTimes[k] = time1 + pathTravelTime + time2
		end
	end
	
	coverBoundMode = CoverBoundMode()
	coverBoundMode.stationsToNodesTimes = stationsToNodesTimes
	coverBoundMode.stationsToPointsTimes = stationsToPointsTimes
	coverBoundMode.pointsToHospitalTimes = pointsToHospitalTimes
	coverBoundMode.nodeSets = nodeSets
	coverBoundMode.nodeSetsStationList = nodeSetsStationList
	coverBoundMode.responseTravelMode = responseTravelMode
	
	return coverBoundMode
end

# mutates: coverBoundMode.distrs, coverBoundMode.transportProb
function setCoverBoundModeDistrs!(coverBoundMode::CoverBoundMode, genConfig::GenConfig)
	distrs = coverBoundMode.distrs # shorthand
	distrs["onSceneDuration"] = genConfig.onSceneDurationDistrRng.d
	distrs["transport"] = genConfig.transportDistrRng.d
	distrs["handoverDuration"] = genConfig.handoverDurationDistrRng.d
	coverBoundMode.transportProb = distrs["transport"].p # probability from bernoulli distribution
end

# mutates: coverBoundMode.distrs, coverBoundMode.transportProb
function setCoverBoundModeDistrs!(coverBoundMode::CoverBoundMode, sim::Simulation)
	genConfig = readGenConfig(sim.inputFiles["callGenConfig"].path)
	setCoverBoundModeDistrs!(coverBoundMode, genConfig)
end

# Create a discrete, non-parametric distribution that approximates the given distribution at the given x values.
# Set lowerBound = true (default) to create distribution that underestimates values, lowerBound = false to overestimate.
# Set removeTrailingZeros = true (default is false) to remove zeros at start and end of cdf of distribution.
function approxDistr(d::Distribution, x::Vector{Float}; lowerBound::Bool = true, removeTrailingZeros::Bool = false)
	x = issorted(x) ? x : sort(x)
	P = [cdf(d,v) for v in x] # probability of sampling value <= v for each v in x
	p = lowerBound ? (vcat(P[2:end],1) - P) : (P - vcat(0,P[1:end-1])) # convert cdf to pdf
	if removeTrailingZeros
		i = 1
		while p[i] == 0 i += 1 end
		j = length(p)
		while p[j] == 0 j-= 1 end
		x = x[i:j]
		p = p[i:j]
	end
	return DiscreteNonParametric(x,p)
end

# Convolute two discrete, non-parametric distributions.
# Requires that x values of the distributions have same step sizes.
function convolute(d1::DiscreteNonParametric, d2::DiscreteNonParametric)
	# shorthand
	x1, p1 = d1.support, d1.p
	x2, p2 = d2.support, d2.p
	n1, n2 = length(x1), length(x2)
	
	if n1 == 1 || n2 == 1
		x = x1 .+ x2
		p = p1 .* p2
		return DiscreteNonParametric(x,p)
	end
	
	# check that distribution of x values in the distributions match
	function getDistributionStepSize(d::DiscreteNonParametric)
		x = d.support
		@assert(issorted(x))
		n = length(x)
		dx = (x[end] - x[1]) / (n - 1)
		@assert(dx > 0)
		for i = 1:n @assert(isapprox(x[i], x[1] + dx*(i-1))) end
		return dx
	end
	dx1 = getDistributionStepSize(d1)
	dx2 = getDistributionStepSize(d2)
	@assert(isapprox(dx1, dx2))
	dx = (dx1 + dx2) / 2
	
	# calculate x for convolution
	xMin = x1[1] + x2[1]
	xMax = x1[end] + x2[end]
	nx = round(Int, (xMax - xMin) / dx) + 1
	x = collect(range(xMin; stop = xMax, length = nx))
	@assert(isapprox((x[end] - x[1]) / (nx - 1), dx))
	
	# calculate p for convolution
	p = zeros(Float, nx) # note: calculating p from p1 * p2' is inefficient
	for i = 1:n1, j = 1:n2
		p[i+j-1] += p1[i] * p2[j]
	end
	@assert(isapprox(sum(p), 1))
	
	return DiscreteNonParametric(x,p)
end

# Calculate distributions needed for lower bound on amb busy duration distribution.
# Mutates: coverBoundMode.approxDistrs
function calcApproxDistrs!(coverBoundMode::CoverBoundMode; pMax::Float = 0.999, dt::Float = 1/60/60/24)
	@assert(0 <= pMax < 1)
	@assert(dt > 0)
	
	function approxDistr2(d::Distribution, pMax::Float, dx::Float; kwargs...)
		@assert(0 <= pMax < 1)
		@assert(dx > 0)
		xMax = quantile(d, pMax)
		xMax = ceil(xMax / dx) * dx
		x = collect(0:dx:xMax)
		return approxDistr(d, x; kwargs...)
	end
	
	kwargs = (lowerBound = true, removeTrailingZeros = true)
	approxDistrs = coverBoundMode.approxDistrs # shorthand
	for (name, distr) in coverBoundMode.distrs
		if name == "transport" continue end
		approxDistrs[name] = approxDistr2(distr, pMax, dt; kwargs...)
	end
	
	approxDistrs["d1"] = deepcopy(approxDistrs["onSceneDuration"]) # approximation of on-scene duration
	approxDistrs["d2"] = convolute(approxDistrs["d1"], approxDistrs["handoverDuration"]) # lower bound approximation of d1 and handover duration
	# Could create a weighted sum of distributions d1 and d2 so that there is only one amb busy duration distribution to sample from,
	# but that is more work (for me) than just having the two distributions and accounting for the relative probability of using each.
	
	return nothing
end

function makeCdf(d::DiscreteNonParametric)
	@assert(issorted(d.support))
	y = copy(d.p)
	for i = 2:length(y) y[i] += y[i-1] end
	@assert(isapprox(y[end], 1))
	return y
end

# Return index i such that xs[i] <= x < xs[i+1], for sorted xs
function binarySearch(xs::Vector{Float}, x::Float)
	# @assert(issorted(xs)) # slow
	n = length(xs)
	
	if x < xs[1] return 0
	elseif x >= xs[end] return n
	end
	
	i = 1; j = n
	while i < j-1
		k = div(i+j, 2)
		xs[k] <= x ? i = k : j = k
	end
	
	return i
end

function binarySearch(d::DiscreteNonParametric, cdf::Vector{Float}, x::Float; xMin::Float = 0.0)
	@assert(length(d.p) == length(cdf))
	@assert(d.support[1] >= xMin) # should this be > instead of >= ?
	i = binarySearch(d.support, x)
	if i == 0 return xMin end
	@assert(1 <= i <= length(d.p))
	return cdf[i]
end

function calcNodeSetsAmbBusyDurationProbs(sim::Simulation, coverBound::CoverBound, targetAmbBusyDuration::Float)
	@assert(targetAmbBusyDuration >= 0)
	
	# shorthand:
	@unpack stations, numStations = sim
	points = sim.demandCoverage.points
	numPoints = length(points)
	nodesPoints = sim.demandCoverage.nodesPoints
	rastersPointDemands = sim.demandCoverage.rastersPointDemands
	
	@assert(sim.demand.numSets == 1) # otherwise would have to calculate cover bound for each demand set
	totalDemand = sum(mode -> mode.arrivalRate, sim.demand.modes) # assumes sim.demand.numSets == 1
	
	stationsNodeSetsP = Dict{Vector{Int}, Vector{Float}}() # stationsNodeSetsP[i][j] = fraction of total demand that can be served by an ambulance that is busy for duration <= targetAmbBusyDuration, by assigning station j to serve node set that has stationList == i
	for coverBoundMode in coverBound.modes, stationList in coverBoundMode.nodeSetsStationList
		get!(stationsNodeSetsP, stationList, zeros(Float, sim.numStations))
	end
	
	for coverBoundMode in coverBound.modes
		# shorthand
		@unpack stationsToPointsTimes, pointsToHospitalTimes, nodeSets, nodeSetsStationList, transportProb = coverBoundMode
		d1 = coverBoundMode.approxDistrs["d1"]
		d2 = coverBoundMode.approxDistrs["d2"]
		numNodeSets = length(nodeSets)
		
		# create cdf for quick look-up of probability values from approximate distributions
		cdf1 = makeCdf(d1)
		cdf2 = makeCdf(d2)
		
		# calculate probabilities for amb being busy for <= targetAmbBusyDuration, for each pair of station and point
		stationsPointsP = zeros(Float, numStations, numPoints) # stationsPointsP[i,j] = probability of amb busy duration (to serve point j) <= targetAmbBusyDuration from station i
		for i = 1:numStations, j = 1:numPoints
			# need to calculate probability that amb will be busy for <= targetAmbBusyDuration while serving point, based on two possible outcomes:
			# - no transport to hospital, look at probability of targetAmbBusyDuration >= durations: response travel + on-scene.
			# - transport to hospital, look at probability of targetAmbBusyDuration >= durations: response travel + on-scene + transport + at-hospital.
			# calculate the above probabilities and multiply by the probability of that outcome, sum for total probability of amb being busy for time <= targetAmbBusyDuration.
			t1 = stationsToPointsTimes[i,j]
			t2 = t1 + pointsToHospitalTimes[j]
			p1 = binarySearch(d1, cdf1, targetAmbBusyDuration - t1) # probability of amb busy duration <= targetAmbBusyDuration, assuming no transport to hospital
			p2 = binarySearch(d2, cdf2, targetAmbBusyDuration - t2) # probability of amb busy duration <= targetAmbBusyDuration, assuming transport to hospital
			stationsPointsP[i,j] = (1-transportProb) * p1 + (transportProb) * p2
		end
		
		# get all demand modes which give the right travel mode (based on response) for this coverBoundMode
		demandModes = DemandMode[]
		for (demandPriority, travelPriority) in sim.responseTravelPriorities
			travelMode = getTravelMode!(sim.travel, travelPriority, sim.startTime) # assumes travel.numSets == 1
			if travelMode == coverBoundMode.responseTravelMode
				demandMode = getDemandMode!(sim.demand, demandPriority, sim.startTime) # assumes demand.numSets == 1
				push!(demandModes, demandMode)
			end
		end
		
		pointDemandFracs = sum(mode -> rastersPointDemands[mode.rasterIndex] * mode.rasterMultiplier, demandModes) / totalDemand # pointDemandFracs[i] is fraction of demand at point i
		
		# group demand points
		for (i,nodeSet) in enumerate(nodeSets)
			stationList = nodeSetsStationList[i]
			nodeSetsP = stationsNodeSetsP[stationList]
			for j in nodeSet # j = node index
				for k in nodesPoints[j] # indices of points for which node j is the nearest node
					for s = 1:numStations
						nodeSetsP[s] += stationsPointsP[s,k] * pointDemandFracs[k]
					end
				end
			end
		end
	end
	
	return stationsNodeSetsP
end

# Get an upper bound on amb busy duration distribution for each
# combination of amb busy duration and number of free ambulances.
function calcAmbBusyDurationProbUpperBounds(sim, coverBound::CoverBound;
		ambBusyDurationsToSample::Vector{Float} = coverBound.ambBusyDurationsToSample,
		numAmbsList::Vector{Int} = collect(1:sim.numStations), pMedianRelax::Bool = true, doPrint::Bool = false)
	
	@assert(issorted(ambBusyDurationsToSample))
	@assert(issorted(numAmbsList))
	
	pMedianOptions = copy(pMedianDefaultOptions)
	pMedianOptions[:x_bin] = !pMedianRelax
	
	numTargets = length(ambBusyDurationsToSample) # shorthand
	n = length(numAmbsList) # shorthand
	ambBusyDurationProbUpperBounds = zeros(Float, numTargets, length(numAmbsList)) # ambBusyDurationProbUpperBounds[i,j] = upper bound on probability of dispatched amb being busy for duration <= ambBusyDurationsToSample[i] for j free ambulances
	for (i,t) in enumerate(ambBusyDurationsToSample)
		doPrint && print("\rTarget amb busy duration $i of $numTargets, 0 ambulances of $n ")
		
		stationsNodeSetsP = calcNodeSetsAmbBusyDurationProbs(sim, coverBound, t)
		stationsP = hcat(values(stationsNodeSetsP)...)
		
		for (j,a) in enumerate(numAmbsList)
			doPrint && print("\rTarget amb busy duration $i of $numTargets, $j ambulances of $n ")
			
			results = Dict()
			solvePMedian(a, -stationsP; options = pMedianOptions, results = results)
			ambBusyDurationProbUpperBounds[i,j] = -results[:cost]
			
			# # Testing: group node sets further, as we can assume that each station will have at least its a^th-to-last preferred station selected.
			# # Maybe only do this if pMedianRelax = false? Otherwise could lead to weaker upper bounds.
			# m = sim.numStations - a + 1
			# stationsReducedNodeSetsP = Dict{Vector{Int}, Vector{Float}}()
			# for (stationList, nodeSetsP) in stationsNodeSetsP
				# x = get!(stationsReducedNodeSetsP, stationList[1:m], zeros(Float, sim.numStations))
				# x .+= nodeSetsP
			# end
			# stationsReducedP = hcat(values(stationsReducedNodeSetsP)...)
			# results = Dict()
			# solvePMedian(a, -stationsReducedP; options = pMedianOptions, results = results)
			# ambBusyDurationProbUpperBounds[i,j] = -results[:cost]
		end
	end
	doPrint && println()
	
	# fix any values in ambBusyDurationProbUpperBounds that are incorrectly larger than they should be (not sure why this happens),
	# can have cases where decreasing the number of ambulances or target amb busy duration can (incorrectly) slightly increase the value from p-median problem.
	p = ambBusyDurationProbUpperBounds # shorthand
	# should have p[i,j] <= p[i+1,j] and p[i,j] <= p[i,j+1]
	m, n = size(p)
	for i = m:-1:1, j = n:-1:1
		if i < m p[i,j] = min(p[i,j], p[i+1,j]) end
		if j < n p[i,j] = min(p[i,j], p[i,j+1]) end
	end
	
	return ambBusyDurationProbUpperBounds
end

# create amb busy duration distribution
# requires fields of coverBound to already be populated: ambBusyDurationsToSample, ambBusyDurationProbUpperBounds
function calcAmbBusyDurationLowerBoundDistrs!(coverBound::CoverBound)
	coverBound.ambBusyDurationLowerBoundDistrs = Sampleable[]
	ambBusyDurations = vcat(0, coverBound.ambBusyDurationsToSample)
	for j = 1:size(coverBound.ambBusyDurationProbUpperBounds,2)
		p1 = Float.(vcat(0, coverBound.ambBusyDurationProbUpperBounds[:,j], 1))
		p2 = p1[2:end] - p1[1:end-1]
		distr = DiscreteNonParametric(ambBusyDurations, p2)
		push!(coverBound.ambBusyDurationLowerBoundDistrs, sampler(distr))
	end
end

# solve MCLP to determine maximum coverage achievable with different numbers of free ambulances (1:sim.numStations)
function calcNumAmbsMaxCoverageFrac(sim::Simulation)
	@assert(sim.demand.numSets == 1)
	@assert(sim.travel.numSets == 1)
	
	# solve MCLP for each possible number of ambulances
	numAmbsMaxCoverageFrac = zeros(Float, sim.numStations) # numAmbsMaxCoverageFrac[i] gives maximum fraction of demand coverage (calculated from the Maximal Coverage Location Problem) achievable with i free ambulances
	# note that having more ambs than stations has no benefit for MCLP; once every station has an amb then max coverage has been achieved
	arrivalRate = sum([mode.arrivalRate for mode in sim.demand.modes]) # assumes demand.numSets == 1
	results = Dict() # for MCLP results
	for i = 1:sim.numStations
		solveMclp!(sim; numAmbs = i, results = results)
		numAmbsMaxCoverageFrac[i] = results[:objVal] / arrivalRate
	end
	
	return numAmbsMaxCoverageFrac
end

function calcQueuedDurationsMaxCoverageFrac(sim::Simulation, coverBound::CoverBound;
	queuedDurationsToSample::Vector{Float} = coverBound.queuedDurationsToSample,
	demandWeights::Dict{Priority,Float} = Dict([p => 1.0 for p in priorities]),
	doPrint::Bool = false)
	
	@assert(queuedDurationsToSample != [])
	@assert(issorted(queuedDurationsToSample))
	@assert(queuedDurationsToSample[1] >= 0)
	
	@assert(all(x -> x >= 0, values(demandWeights)))
	@assert(sum(values(demandWeights)) > 0)
	for p in keys(demandWeights)
		if demandWeights[p] == 0
			delete!(demandWeights, p)
		end
	end
	demandPriorities = collect(keys(demandWeights))
	
	@assert(sim.demand.numSets == 1)
	@assert(sim.travel.numSets == 1)
	
	numStations = sim.numStations
	numPoints = length(sim.demandCoverage.points)
	arrivalRate = sum(mode -> mode.arrivalRate, sim.demand.modes) # assumes demand.numSets == 1
	
	# get data for each demand priority
	# cannot necessarily combine these (e.g. combine by common travel modes), as different demand priorities may have different target response times
	demandPrioritiesData = Dict()
	currentTime = sim.startTime
	for demandPriority in demandPriorities
		responseTravelPriority = sim.responseTravelPriorities[demandPriority]
		demandMode = getDemandMode!(sim.demand, demandPriority, currentTime) # assumes demand.numSets == 1
		travelMode = getTravelMode!(sim.travel, responseTravelPriority, sim.startTime) # assumes sim.travel.numSets == 1
		demandPriorityData = Dict()
		demandPriorityData[:coverBoundMode] = coverBound.modes[coverBound.modeLookup[travelMode.index]]
		demandPriorityData[:rasterPointDemands] = sim.demandCoverage.rastersPointDemands[demandMode.rasterIndex] * demandMode.rasterMultiplier * demandWeights[demandPriority]
		demandPrioritiesData[demandPriority] = demandPriorityData
	end
	
	maxCoverageFracs = zeros(Float, length(queuedDurationsToSample)) # maxCoverageFracs[i] = mclp value for queuedDurationsToSample[i]
	# stationsNumAmbsList = []
	for (i,queuedDuration) in enumerate(queuedDurationsToSample)
		doPrint && print("\rqueuedDuration $i of $(length(queuedDurationsToSample))")
		
		# calculate demand covered by different station sets
		pointData = Dict{Vector{Bool}, Float}() # pointData[stationsCoverPoint] = pointDemand
		for demandPriority in demandPriorities
			coverTime = sim.targetResponseDurations[Int(demandPriority)] - (coverBound.dispatchDelay + queuedDuration)
			if coverTime > 0
				demandPriorityData = demandPrioritiesData[demandPriority]
				rasterPointDemands = demandPriorityData[:rasterPointDemands]
				stationsCoverPoints = demandPriorityData[:coverBoundMode].stationsToPointsTimes .<= coverTime # would need to change this if dispatch delay is not deterministic
				for j = 1:numPoints
					stationsCoverPoint = stationsCoverPoints[:,j] # stationSet[i] = true if station i covers point j
					get!(pointData, stationsCoverPoint, 0.0)
					pointData[stationsCoverPoint] += rasterPointDemands[j]
				end
			end
		end
		pointStations = findall.(collect(keys(pointData)))
		pointDemands = collect(values(pointData))
		# point j has demand pointDemands[j] and is covered by stations pointStations[j]
		
		if isempty(pointStations) || (length(pointStations) == 1 && !any(pointStations[1]))
			break # cannot reach any points in time
		end
		
		# solve mclp
		results = Dict()
		stationsNumAmbs = solveMclp(1, pointDemands, pointStations, results = results)
		maxCoverageFracs[i] = results[:objVal] / arrivalRate
		# push!(stationsNumAmbsList, stationsNumAmbs)
	end
	
	# # if call was queued for no time, coverage should match that of MCLP for one amb and full target response time (minus dispatch delay)
	# if queuedDurationsToSample[1] == 0 @assert(isapprox(maxCoverageFracs[1], coverBound.numAmbsMaxCoverageFrac[1])) end
	
	return maxCoverageFracs # max probability of covering call that has been queued for queuedDurationsToSample[i] is maxCoverageFracs[i]
end

function initCoverBoundSim(; numReps::Int = 1, numAmbs::Int = 0, warmUpDuration::Float = 0.0, minLastCallArrivalTime::Float = nullTime,
		interarrivalTimeDistrRng::Union{DistrRng,Nothing} = nothing,
		arrivalRate::Float = 0.0, interarrivalTimeSeed::Int = 0, # use these if interarrivalTimeDistrRng == nothing
		ambBusyDurationSeed::Int = 1)
	
	@assert(numReps >= 1)
	@assert(numAmbs >= 1)
	@assert(warmUpDuration >= 0)
	@assert(minLastCallArrivalTime > warmUpDuration)
	
	if interarrivalTimeDistrRng == nothing
		@assert(arrivalRate > 0)
		interarrivalTimeDistrRng = DistrRng(Exponential(1/arrivalRate), seed = interarrivalTimeSeed)
	else
		@assert(isa(interarrivalTimeDistrRng.d, Exponential))
		@assert(mean(interarrivalTimeDistrRng.d) != Inf)
	end
	
	ambBusyDurationRng = MersenneTwister(ambBusyDurationSeed)
	
	coverBoundSim = CoverBoundSim()
	@pack! coverBoundSim = numReps, numAmbs, warmUpDuration, minLastCallArrivalTime, interarrivalTimeDistrRng, ambBusyDurationRng
	
	return coverBoundSim
end

# mutates: coverBoundSim.interarrivalTimeDistrRng, coverBoundSim.ambBusyDurationRng
function simulateRepCoverBound!(coverBound::CoverBound)
	cbs = coverBound.sim # shorthand
	@assert(cbs.warmUpDuration > 0)
	@assert(cbs.interarrivalTimeDistrRng != nothing)
	@assert(cbs.numAmbs >= 1)
	@assert(0 <= cbs.warmUpDuration < cbs.minLastCallArrivalTime)
	
	numStations = length(coverBound.ambBusyDurationLowerBoundDistrs)
	@assert(numStations >= 1)
	
	# generate call arrival times
	# could do this within sim loop, but will do here for simplicity
	callArrivalTimes = Float[]
	t = 0.0
	while t < cbs.minLastCallArrivalTime
		t += rand(cbs.interarrivalTimeDistrRng)
		push!(callArrivalTimes, t)
	end
	numCalls = length(callArrivalTimes)
	@assert(numCalls >= 1)
	
	# run simple sim, count number of free ambulances at moment of each dispatch
	currentTime = 0.0
	eventList = Event[]
	addEvent!(eventList, form = callArrives, time = callArrivalTimes[1])
	numFreeAmbs = cbs.numAmbs
	numFreeAmbsCount = zeros(Int, cbs.numAmbs) # numFreeAmbsCount[i] = how many times calls arrive when there are i free ambulances
	numQueuedCalls = 0
	earliestQueuedCallIndex = nullIndex # of the currently queued calls, index of the call that arrived earliest
	queuedCallDurations = Float[] # duration for which calls (past and present) were queued
	callIndex = 1
	while !isempty(eventList)
		event = getNextEvent!(eventList)
		currentTime = event.time
		if event.form == callArrives
			event = addEvent!(eventList, form = considerDispatch, time = callArrivalTimes[callIndex] + coverBound.dispatchDelay)
			event.callIndex = callIndex
			callIndex += 1
			if callIndex <= numCalls
				event = addEvent!(eventList, form = callArrives, time = callArrivalTimes[callIndex])
			end
			
		elseif event.form == considerDispatch
			if numFreeAmbs == 0
				numQueuedCalls += 1
				if earliestQueuedCallIndex == nullIndex
					earliestQueuedCallIndex = event.callIndex
				end
			else
				addEvent!(eventList; form = ambDispatched, time = currentTime)
			end
			
		elseif event.form == ambDispatched
			recordStats = event.time >= cbs.warmUpDuration
			if earliestQueuedCallIndex == nullIndex # no queued calls
				recordStats && (numFreeAmbsCount[numFreeAmbs] += 1)
			else
				# assume that queued calls are responded to in first-in first-out order
				queuedDuration = event.time - (callArrivalTimes[earliestQueuedCallIndex] + coverBound.dispatchDelay)
				@assert(queuedDuration >= 0)
				recordStats && push!(queuedCallDurations, queuedDuration)
				earliestQueuedCallIndex = numQueuedCalls == 0 ? nullIndex : earliestQueuedCallIndex + 1
			end
			
			# generate amb busy duration
			t = rand(cbs.ambBusyDurationRng, coverBound.ambBusyDurationLowerBoundDistrs[min(numFreeAmbs, numStations)])
			addEvent!(eventList; form = ambBecomesFree, time = currentTime + t)
			
			numFreeAmbs -= 1
		elseif event.form == ambBecomesFree
			numFreeAmbs += 1
			if numQueuedCalls >= 1
				addEvent!(eventList; form = ambDispatched, time = currentTime)
				numQueuedCalls -= 1
			end
		else error("Unrecognised event.")
		end
	end
	
	return numCalls, numFreeAmbsCount, queuedCallDurations # numCalls = sum(numFreeAmbsCount) + length(queuedCallDurations)
end

# mutates: coverBound.sim.reps, coverBound.sim.bound
function simulateCoverBound!(coverBound::CoverBound)
	coverBound.sim.reps = [CoverBoundSimRep() for i = 1:coverBound.sim.numReps]
	for rep in coverBound.sim.reps
		rep.numCalls, rep.numFreeAmbsCount, rep.queuedCallDurations = simulateRepCoverBound!(coverBound)
	end
	calcCoverBound!(coverBound)
	return coverBound.sim.bound
end

# mutates: coverBound.sim.bound, rep.bound for rep in coverBound.sim.reps
function calcCoverBound!(coverBound::CoverBound)
	cb = coverBound # shorthand
	for (i,rep) in enumerate(cb.sim.reps)
		@unpack numCalls, numFreeAmbsCount, queuedCallDurations = rep
		n = min(length(cb.numAmbsMaxCoverageFrac), length(numFreeAmbsCount))
		coverageFracTotal = sum(cb.numAmbsMaxCoverageFrac[1:n] .* numFreeAmbsCount[1:n]) + cb.numAmbsMaxCoverageFrac[end] * sum(numFreeAmbsCount[n+1:end])
		if cb.accountForQueuedDurations
			# account for duration that calls were queued
			for queuedCallDuration in queuedCallDurations
				j = binarySearch(cb.queuedDurationsToSample, queuedCallDuration)
				if j == 0 # queuedCallDuration < queuedDurationsToSample[1], so need to assume call was queued for no time
					coverageFracTotal += cb.numAmbsMaxCoverageFrac[1] # maximum coverage from one ambulance
				else # j > 0
					coverageFracTotal += cb.queuedDurationsMaxCoverageFrac[j] # coverage accounting for queued time
				end
			end
		else
			# original cover bound, ignores duration for which calls were queued and gives these calls the maximum coverage from one ambulance
			coverageFracTotal += cb.numAmbsMaxCoverageFrac[1] * length(queuedCallDurations)
		end
		rep.bound = coverageFracTotal / numCalls
	end
	
	bounds = [rep.bound for rep in cb.sim.reps]
	cb.sim.bound = MeanAndHalfWidth(mean(bounds), tDistrHalfWidth(bounds))
	
	return cb.sim.bound
end
