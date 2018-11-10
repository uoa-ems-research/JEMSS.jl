# Maximum Expected Coverage Location Problem (MEXCLP)

"""
	solveMexclp!(sim::Simulation;
		numAmbs::Int = length(sim.ambulances),
		busyFraction::Float = 0.5,
		demandWeights::Dict{Priority,Float} = Dict([p => 1.0 for p in instances(Priority)]))
Solves the Maximum Expected Coverage Location Problem (MEXCLP) for `sim` and returns the number of ambulances to assign to each station, also the converse is returned - a station index for each ambulance.
The problem assumes that all ambulances are equivalent.
Requires data for demand and demand coverage to already be set in `sim`; see functions [`initDemand!`](@ref), [`initDemandCoverage!`](@ref).

# Keyword arguments
- `numAmbs` is the number of ambulances to solve for, must be >= 1
- `busyFraction` is the fraction of time that ambulances are busy; should be within [0,1] though this is not enforced
- `demandWeights` is the weight to apply to each demand priority for the objective function
"""
function solveMexclp!(sim::Simulation;
	numAmbs::Int = length(sim.ambulances),
	busyFraction::Float = 0.5,
	demandWeights::Dict{Priority,Float} = Dict([p => 1.0 for p in instances(Priority)]))
	
	@assert(numAmbs >= 1, "need at least 1 ambulance for mexclp")
	@assert(sim.travel.numSets == 1) # otherwise need to solve mexclp for each travel set?
	@assert(sim.demand.numSets != nullIndex, "no demand data, cannot solve mexclp")
	@assert(sim.demand.numSets == 1) # otherwise need to solve mexclp for each demand set?
	
	# check that demand coverage data is set
	@assert(!isempty(sim.demandCoverTimes))
	@assert(!isempty(sim.demandCoverage.points))
	
	# shorthand
	numStations = length(sim.stations)
	currentTime = sim.startTime
	
	# get demand point coverage data
	pointData = Dict{Vector{Int},Float}() # pointData[pointStations] = pointDemand
	demandPriorities = setdiff([instances(Priority)...], [nullPriority])
	for demandPriority in demandPriorities
		if !haskey(demandWeights, demandPriority) || demandWeights[demandPriority] == 0
			continue
		end
		
		pointsCoverageMode = getPointsCoverageMode!(sim, demandPriority, currentTime)
		demandMode = getDemandMode!(sim.demand, demandPriority, currentTime)
		pointSetsDemands = getPointSetsDemands!(sim, demandPriority, currentTime; pointsCoverageMode = pointsCoverageMode) * demandMode.rasterMultiplier * demandWeights[demandPriority]
		
		for (i,stationSet) in enumerate(pointsCoverageMode.stationSets)
			get!(pointData, stationSet, 0.0)
			pointData[stationSet] += pointSetsDemands[i]
		end
	end
	pointStations = collect(keys(pointData))
	pointDemands = collect(values(pointData))
	# point j has demand pointDemands[j] and is covered by stations pointStations[j]
	numPoints = length(pointDemands) # shorthand
	
	# calculate benefit of covering each point with a kth ambulance
	marginalBenefit = (busyFraction.^[0:numAmbs-1;])*(1-busyFraction) # point cover benefit values, for single demand
	pointCoverCountBenefit = [pointDemands[j] * marginalBenefit for j = 1:numPoints] # pointCoverCountBenefit[j][k] is the benefit of adding a kth ambulance to cover point j
	
	#######
	# IP
	
	# shorthand:
	a = numAmbs
	s = numStations
	p = numPoints
	
	model = Model()
	
	if !all(v -> issorted(v, rev=true), pointCoverCountBenefit)
		setsolver(model, GLPKSolverMIP(presolve=true)) # GLPKSolverMIP solves faster than CbcSolver for this formulation
		
		@variables(model, begin
			(x[i=1:s] >= 0, Int) # x[i] = number of ambulances assigned to station i
			(y[j=1:p,k=1:a], Bin) # y[j,k] = true if point j is covered by at least k ambulances
		end)
		
		@constraints(model, begin
			(useAllAmbs, sum(x) == a)
			(pointCoverCount[j=1:p], sum(y[j,:]) == sum(x[pointStations[j]]))
			(pointCoverOrder[j=1:p, k=1:(a-1)], y[j,k] >= y[j,k+1])
		end)
	else
		# pointCoverCountBenefit[j] is non-increasing (because busyFraction >= 0), so:
		# - can have y variables as non-binary (but still constrained to be between 0 and 1)
		# - can leave out the constraint 'pointCoverOrder'
		# - could have x variables as non-integer, though if there are multiple solutions then x values might not be naturally integer, so will leave the integer constraint
		
		setsolver(model, CbcSolver()) # CbcSolver solves faster than GLPKSolverMIP for this formulation
		# using CbcSolver here would be faster with presolve, but using this causes a line print when solving
		
		@variables(model, begin
			(x[i=1:s] >= 0, Int) # x[i] = number of ambulances assigned to station i
			(0 <= y[j=1:p,k=1:a] <= 1) # y[j,k] = true if point j is covered by at least k ambulances
		end)
		
		@constraints(model, begin
			(useAllAmbs, sum(x) == a)
			(pointCoverCount[j=1:p], sum(y[j,:]) == sum(x[pointStations[j]]))
		end)
	end
	
	@expressions(model, begin
		# (expectedPointCoverage[j=1:p], sum(y[j,k] * pointCoverCountBenefit[j][k] for k=1:a))
		(expectedCoverage, sum(y[j,k] * pointCoverCountBenefit[j][k] for j=1:p, k=1:a))
	end)
	
	@objective(model, :Max, expectedCoverage)
	
	solve(model)
	
	# if IP is slow, try:
	# - if pointCoverCountBenefit[j] is decreasing, then can change y variables to be non-binary (but still constrained to be between 0 and 1), and remove the constraint 'pointCoverOrder'; can also remove constraint that x is integer (assuming that there is one unique solution...)
	
	# extract solution
	stationsNumAmbs = convert(Vector{Int}, round.(getvalue(x))) # solution; stationsNumAmbs[i] gives number of ambulances to allocate to station i
	pointSlotCover = convert(Array{Bool,2}, round.(getvalue(y))) # pointSlotCover[j,k] = true if point j is covered by at least k ambulances
	
	if checkMode
		# check constraints
		@assert(sum(stationsNumAmbs) == a) # useAllAmbs constraint
		@assert(all(j -> sum(pointSlotCover[j,:]) == sum(stationsNumAmbs[pointStations[j]]), 1:p)) # pointCoverCount constraint
		@assert(all(j -> issorted(pointSlotCover[j], rev=true), 1:p)) # pointCoverOrder constraint
	end
	
	# convert to a deployment policy
	depol = Depol() # depol[i] gives the station index for ambulance i
	for (stationIndex, numAmbs) in enumerate(stationsNumAmbs), i = 1:numAmbs
		push!(depol, stationIndex)
	end
	
	return stationsNumAmbs, depol
end
