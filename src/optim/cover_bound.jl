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

function calcCoverBound!(sim::Simulation;
		coverBound::Union{CoverBound,Nothing} = nothing,
		coverBoundSim::Union{CoverBoundSim,Nothing} = nothing, # need this if coverBound.sim is not set
		serviceDurationsToSample::Vector{Float} = [],
		doPrint::Bool = false)
	
	@assert(serviceDurationsToSample != [] && all(serviceDurationsToSample .> 0))
	@assert(issorted(serviceDurationsToSample))
	
	if coverBound == nothing
		coverBound = initCoverBound(sim, doPrint = doPrint)
		@assert(coverBoundSim != nothing)
		coverBound.sim = coverBoundSim
	end
	cb = coverBound # shorthand
	cbs = coverBoundSim = coverBound.sim # shorthand
	
	if cb.serviceDurationsToSample != serviceDurationsToSample || length(cb.serviceDurationsToSample) != size(cb.serviceProbUpperBounds, 1)
		cb.serviceDurationsToSample = serviceDurationsToSample
		cb.serviceProbUpperBounds = calcServiceProbUpperBounds(sim, cb, doPrint = doPrint)
	else
		# user should ensure that serviceProbUpperBounds are for serviceDurationsToSample
		doPrint && println("Using existing coverBound.serviceProbUpperBounds.")
	end
	
	calcServiceDurationLowerBoundDistrs!(cb) # create service duration distributions
	
	# set length of numAmbsMaxCoverageFrac (set in initCoverBound function) to match maximum(numAmbs, numStations)
	x = cb.numAmbsMaxCoverageFrac # shorthand
	n = cbs.numAmbs # shorthand
	numAmbsMaxCoverageFrac = n <= length(x) ? x[1:n] : vcat(x, fill(x[end], n - length(x)))
	
	# simulate
	repsNumFreeAmbsFrac = simulateCoverBound!(cb)
	v = [sum(numAmbsMaxCoverageFrac .* fracs) for fracs in repsNumFreeAmbsFrac]
	cbs.bound = MeanAndHalfWidth(mean(v), tDistrHalfWidth(v)) # cover bound result
	
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
	distrs["dispatchDelay"] = genConfig.dispatchDelayDistrRng.d
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

# Calculate distributions needed for lower bound on service time distrubution.
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
	
	approxDistrs["d1"] = convolute(approxDistrs["dispatchDelay"], approxDistrs["onSceneDuration"]) # lower bound approximation of dispatch delay and on-scene duration
	approxDistrs["d2"] = convolute(approxDistrs["d1"], approxDistrs["handoverDuration"]) # lower bound approximation of d1 and handover duration
	# Could create a weighted sum of distributions d1 and d2 so that there is only one service duration distribution to sample from,
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

function calcNodeSetsServiceProbs(sim::Simulation, coverBound::CoverBound, targetServiceDuration::Float)
	@assert(targetServiceDuration >= 0)
	
	# shorthand:
	@unpack stations, numStations = sim
	points = sim.demandCoverage.points
	numPoints = length(points)
	nodesPoints = sim.demandCoverage.nodesPoints
	rastersPointDemands = sim.demandCoverage.rastersPointDemands
	
	@assert(sim.demand.numSets == 1) # otherwise would have to calculate cover bound for each demand set
	totalDemand = sum(mode -> mode.arrivalRate, sim.demand.modes) # assumes sim.demand.numSets == 1
	
	stationsNodeSetsP = Dict{Vector{Int}, Vector{Float}}() # stationsNodeSetsP[i][j] = fraction of total demand that can be served within targetServiceDuration by assigning station j to serve node set that has stationList == i
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
		
		# calculate probability of station serving demand point in time, for each pair of station and point
		stationsPointsP = zeros(Float, numStations, numPoints) # stationsPointsP[i,j] = probability of serving point j within targetServiceDuration from station i
		for i = 1:numStations, j = 1:numPoints
			# need to calculate probability of serving point within targetServiceDuration based on two possible outcomes:
			# - no transport to hospital, look at probability of targetServiceDuration >= durations: dispatch delay + response travel + on-scene.
			# - transport to hospital, look at probability of targetServiceDuration >= durations: dispatch delay + response travel + on-scene + transport + at-hospital.
			# calculate the above probabilities and multiply by the probability of that outcome, sum for total probability of service within targetServiceDuration.
			t1 = stationsToPointsTimes[i,j]
			t2 = t1 + pointsToHospitalTimes[j]
			p1 = binarySearch(d1, cdf1, targetServiceDuration - t1) # probability of serving point within targetServiceDuration, assuming no transport to hospital
			p2 = binarySearch(d2, cdf2, targetServiceDuration - t2) # probability of serving point within targetServiceDuration, assuming transport to hospital
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

# Get an upper bound on service duration distribution for each
# combination of service duration and number of free ambulances.
function calcServiceProbUpperBounds(sim, coverBound::CoverBound;
		serviceDurationsToSample::Vector{Float} = coverBound.serviceDurationsToSample,
		numAmbsList::Vector{Int} = collect(1:sim.numStations), pMedianRelax::Bool = true, doPrint::Bool = false)
	
	pMedianOptions = copy(pMedianDefaultOptions)
	pMedianOptions[:x_bin] = !pMedianRelax
	
	numTargets = length(serviceDurationsToSample) # shorthand
	n = length(numAmbsList) # shorthand
	serviceProbUpperBounds = zeros(Float, numTargets, length(numAmbsList)) # serviceProbUpperBounds[i,j] = upper bound on probability of serving next call within serviceDurationsToSample[i] for numAmbsList[j] free ambulances
	for (i,t) in enumerate(serviceDurationsToSample)
		doPrint && print("\rTarget service duration $i of $numTargets, 0 ambulances of $n ")
		
		stationsNodeSetsP = calcNodeSetsServiceProbs(sim, coverBound, t)
		stationsP = hcat(values(stationsNodeSetsP)...)
		
		for (j,a) in enumerate(numAmbsList)
			doPrint && print("\rTarget service duration $i of $numTargets, $j ambulances of $n ")
			
			results = solvePMedian(-stationsP, a; options = pMedianOptions)
			serviceProbUpperBounds[i,j] = -results[:cost]
			
			# # Testing: group node sets further, as we can assume that each station will have at least its a^th-to-last preferred station selected.
			# # Maybe only do this if pMedianRelax = false? Otherwise could lead to weaker upper bounds.
			# m = sim.numStations - a + 1
			# stationsReducedNodeSetsP = Dict{Vector{Int}, Vector{Float}}()
			# for (stationList, nodeSetsP) in stationsNodeSetsP
				# x = get!(stationsReducedNodeSetsP, stationList[1:m], zeros(Float, sim.numStations))
				# x .+= nodeSetsP
			# end
			# stationsReducedP = hcat(values(stationsReducedNodeSetsP)...)
			# results = solvePMedian(-stationsReducedP, a; options = pMedianOptions)
			# serviceProbUpperBounds[i,j] = -results[:cost]
		end
	end
	doPrint && println()
	
	return serviceProbUpperBounds
end

# create service duration distribution
# requires fields of coverBound to already be populated: serviceDurationsToSample, serviceProbUpperBounds
function calcServiceDurationLowerBoundDistrs!(coverBound::CoverBound)
	coverBound.serviceDurationLowerBoundDistrs = Sampleable[]
	serviceDurations = vcat(0, coverBound.serviceDurationsToSample)
	for j = 1:size(coverBound.serviceProbUpperBounds,2)
		p1 = Float.(vcat(0, coverBound.serviceProbUpperBounds[:,j], 1))
		p2 = p1[2:end] - p1[1:end-1]
		distr = DiscreteNonParametric(serviceDurations, p2)
		push!(coverBound.serviceDurationLowerBoundDistrs, sampler(distr))
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
		solveMexclp!(sim; numAmbs = i, busyFraction = 0.0, results = results) # MEXCLP with busyFraction = 0 is same as MCLP
		numAmbsMaxCoverageFrac[i] = results[:objVal] / arrivalRate
	end
	
	return numAmbsMaxCoverageFrac
end

function initCoverBoundSim(; numReps::Int = 1, numAmbs::Int = 0, warmUpDuration::Float = 0.0, minLastCallArrivalTime::Float = nullTime,
		interarrivalTimeDistrRng::Union{DistrRng,Nothing} = nothing,
		arrivalRate::Float = 0.0, interarrivalTimeSeed::Int = 0, # use these if interarrivalTimeDistrRng == nothing
		serviceDurationSeed::Int = 1)
	
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
	
	serviceDurationRng = MersenneTwister(serviceDurationSeed)
	
	coverBoundSim = CoverBoundSim()
	@pack! coverBoundSim = numReps, numAmbs, warmUpDuration, minLastCallArrivalTime, interarrivalTimeDistrRng, serviceDurationRng
	
	return coverBoundSim
end

# mutates: coverBoundSim.interarrivalTimeDistrRng, coverBoundSim.serviceDurationRng
function simulateRepCoverBound!(coverBound::CoverBound)
	cbs = coverBound.sim # shorthand
	@assert(cbs.warmUpDuration > 0)
	@assert(cbs.interarrivalTimeDistrRng != nothing)
	@assert(cbs.numAmbs >= 1)
	@assert(0 <= cbs.warmUpDuration < cbs.minLastCallArrivalTime)
	
	numStations = length(coverBound.serviceDurationLowerBoundDistrs)
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
	
	# run simple sim, count number of free ambulances at moment of each call arrival
	currentTime = 0.0
	eventList = Event[]
	addEvent!(eventList, form = callArrives, time = callArrivalTimes[1])
	numFreeAmbs = cbs.numAmbs
	numFreeAmbsCounts = zeros(Int, cbs.numAmbs) # numFreeAmbsCounts[i] = how many times calls arrive when there are i free ambulances
	numQueuedCalls = 0
	callIndex = 1
	while !isempty(eventList)
		event = getNextEvent!(eventList)
		currentTime = event.time
		if event.form == callArrives
			if numFreeAmbs == 0
				numQueuedCalls += 1
			else
				addEvent!(eventList; form = ambDispatched, time = currentTime)
			end
			callIndex += 1
			if callIndex <= numCalls
				event = addEvent!(eventList, form = callArrives, time = callArrivalTimes[callIndex])
			end
		elseif event.form == ambDispatched
			if event.time >= cbs.warmUpDuration
				numFreeAmbsCounts[numFreeAmbs] += 1 # use this if cover bound assumes at least one amb available at call arrival (whether queued or not)
			end
			
			# generate service duration
			t = rand(cbs.serviceDurationRng, coverBound.serviceDurationLowerBoundDistrs[min(numFreeAmbs, numStations)])
			addEvent!(eventList; form = ambBecomesFree, time = currentTime + t)
			
			numFreeAmbs -= 1
		elseif event.form == ambBecomesFree
			numFreeAmbs += 1
			if numQueuedCalls >= 1
				addEvent!(eventList; form = ambDispatched, time = currentTime)
				numQueuedCalls -= 1
				# Note that the cover bound ignores the duration for which calls are queued, effectively acting as if they were queued for no time. Otherwise would need to solve MCLP for different response time targets to account for call queueing duration. The lower bound is still valid but may be weaker than it could be, especially if calls are queued often.
			end
		else error("Unrecognised event.")
		end
	end
	
	numFreeAmbsFrac = numFreeAmbsCounts ./ sum(numFreeAmbsCounts)
	
	return numFreeAmbsFrac
end

# mutates: coverBound.sim.reps
function simulateCoverBound!(coverBound::CoverBound)
	coverBound.sim.reps = [CoverBoundSimRep() for i = 1:coverBound.sim.numReps]
	for rep in coverBound.sim.reps
		rep.numFreeAmbsFrac = simulateRepCoverBound!(coverBound)
	end
	return [rep.numFreeAmbsFrac for rep in coverBound.sim.reps]
end
