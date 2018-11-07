# Maximum Expected Coverage Location Problem (MEXCLP)
# This problem assumes that ambulances are homogeneous.

# Solve the MEXCLP problem for the simulation, return solution as a deployment policy
# Requires travel and demand to be temporally static.
# If sim.demand is not already set, will set from kwarg 'demand', otherwise from reading file 'demandFilename'
function solveMexclp!(sim::Simulation;
	demandCoverTimes = Dict([p => 8/24/60 for p in instances(Priority)]), busyFraction::Float = 0.5,
	demand::Union{Demand,Void} = nothing, demandFilename::String = "",
	rasterCellNumPointRows::Int = 1, rasterCellNumPointCols::Int = 1)
	
	warn("will just consider high priority demand, for now.")
	demandPriority = highPriority
	
	@assert(sim.travel.numSets == 1) # otherwise need to solve mexclp for each travel set?
	
	# shorthand
	numAmbs = length(sim.ambulances)
	numStations = length(sim.stations)
	currentTime = sim.startTime
	
	# set demand
	if sim.demand.numSets == nullIndex
		if demand != nothing
			sim.demand = demand
		elseif demandFilename != ""
			sim.demand = readDemandFile(demandFilename)
		else
			warn("No demand data, cannot solve mexclp.")
			return
		end
	end
	@assert(sim.demand.numSets == 1) # otherwise need to solve mexclp for each demand set?
	
	# initialise demand coverage data
	sim.demandCoverTimes = demandCoverTimes
	initDemandCoverage!(sim; rasterCellNumRows = rasterCellNumPointRows, rasterCellNumCols = rasterCellNumPointCols)
	
	# get demand point coverage data
	pointsCoverageMode = getPointsCoverageMode!(sim, demandPriority, currentTime)
	demandMode = getDemandMode!(sim.demand, demandPriority, currentTime)
	pointSetsDemands = getPointSetsDemands!(sim, demandPriority, currentTime; pointsCoverageMode = pointsCoverageMode) * demandMode.rasterMultiplier
	
	# shorthand:
	points = pointsCoverageMode.pointSets # will refer to 'point sets' as just 'points' from now on
	numPoints = length(points)
	pointDemands = pointSetsDemands
	pointStations = pointsCoverageMode.stationSets
	
	# points[j] has demand pointDemands[j] and is covered by stations pointStations[j]
	
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
	
	setsolver(model, GLPKSolverMIP(presolve=true))
	
	@variables(model, begin
		(x[i=1:s] >= 0, Int) # x[i] = number of ambulances assigned to station i
		(y[j=1:p,k=1:a], Bin) # y[j,k] = true if point j is covered by at least k ambulances
	end)
	
	@constraints(model, begin
		(useAllAmbs, sum(x) == a)
		(pointCoverCount[j=1:p], sum(y[j,:]) == sum(x[pointStations[j]]))
		(pointCoverOrder[j=1:p, k=1:(a-1)], y[j,k] >= y[j,k+1])
	end)
	
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
	depol = Depol()
	for (stationIndex, numAmbs) in enumerate(stationsNumAmbs), i = 1:numAmbs
		push!(depol, stationIndex)
	end
	
	return depol # depol[i] gives the station index for ambulance i
end
