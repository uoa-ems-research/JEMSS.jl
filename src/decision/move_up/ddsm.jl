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

# Dynamic Double Standard Model, with some modifications
# From paper: "A dynamic model and parallel tabu search heuristic for real-time ambulance relocation"

# Modifications (compared with original ddsm):
# - Added slack variables to coverage constraints, to avoid infeasibility. Cost of slack per unit is given by slackWeight.
# - Redeployment cost only depends on travel time, not on: ambulance history, avoiding round trips.

# initialise data relevant to move up
function initDdsm!(sim::Simulation;
	alpha::Float = 0.5, travelTimeCost::Float = 50.0, slackWeight::Float = 1000.0,
	coverTimeDemandPriorities::Vector{Priority} = [highPriority, lowPriority],
	options::Dict{Symbol,Any} = Dict{Symbol,Any}())
	# coverTimeDemandPriorities[i] is demand priority for coverTimes[i], where coverTimes[1] and [2] are the targets for ddsm
	
	debugMode && @warn("Have not finished implementing ddsm.")
	debugMode && @warn("Have allowed redeploying after mission, though original DDSM only redeployed on call arrival.")
	debugMode && @info("Should change variable `alpha` to be more descriptive.")
	debugMode && @info("Should try Gurobi with presolve.") # Presolve=0
	debugMode && @warn("Need to check which combinations of bin and float variables give correct solutions.")
	debugMode && @info("Should test solving with IP gap > 0.")
	debugMode && @info("Should move `options` into config xml.")
	
	@assert(0 <= alpha <= 1)
	@assert(length(coverTimeDemandPriorities) == 2)
	
	# initialise demand and demand coverage data if not already initialised
	sim.demand.initialised || initDemand!(sim)
	sim.demandCoverage.initialised || initDemandCoverage!(sim)
	demand = sim.demand # shorthand
	
	# @show([[demandMode.rasterIndex for demandMode in demand.modes[demand.modeLookup[i,:]]] for i = 1:demand.numSets])
	@assert(all([demandMode.rasterIndex for demandMode in demand.modes[demand.modeLookup[i,:]]] |> unique |> length == 1 for i = 1:demand.numSets)) # distribution of demand must be same for all demand priorities, in any given time period
	
	# two target cover times for ddsm
	coverTimes = [sim.demandCoverage.coverTimes[p] for p in coverTimeDemandPriorities]
	@assert(length(coverTimes) == 2)
	@assert(coverTimes[1] < coverTimes[2]) # change to warning?
	
	# set options
	ddsmd = sim.moveUpData.ddsmData # shorthand
	ddsmd.options[:solver] = "cbc" # can be slower than glpk, but more reliable for some reason
	ddsmd.options[:v] = v"1"
	ddsmd.options[:z_var] = false
	ddsmOptions!(sim, options)
	
	ddsmd.alpha = alpha
	ddsmd.travelTimeCost = travelTimeCost
	ddsmd.slackWeight = slackWeight
	ddsmd.coverTimeDemandPriorities = coverTimeDemandPriorities # coverTimeDemandPriorities[i] is demand priority for coverTimes[i]
	ddsmd.coverTimes = coverTimes
end

# change the options for ddsm
function ddsmOptions!(sim, options::Dict{Symbol,Any})
	currentOptions = sim.moveUpData.ddsmData.options # shorthand
	options = merge!(currentOptions, options) # update currentOptions with options
	
	@assert(in(options[:solver], ["cbc", "glpk", "gurobi"]))
	@assert(typeof(options[:v]) == VersionNumber)
	@assert(typeof(options[:z_var]) == Bool)
	
	if options[:v] == v"1"
		options[:x_bin] = true
		options[:y11_bin] = true
		options[:y12_bin] = true
		options[:y2_bin] = true
	elseif options[:v] == v"2"
		options[:x_bin] = true
		options[:y11_bin] = true
		options[:y12_bin] = false
		options[:y2_bin] = false
	elseif options[:v] == v"3"
		options[:x_bin] = false
		options[:y11_bin] = true
		options[:y12_bin] = false
		options[:y2_bin] = false
	elseif options[:v] == v"4"
		options[:x_bin] = false
		options[:y11_bin] = true
		options[:y12_bin] = true
		options[:y2_bin] = true
	elseif options[:v] == v"5"
		# version to customise which variables are bin
		@assert(all(key -> haskey(options, key), [:x_bin, :y11_bin, :y12_bin, :y2_bin]))
	else
		# only have x as bin, though solution is not guaranteed to be correct
		# might be a good bound?
		options[:x_bin] = true
		options[:y11_bin] = false
		options[:y12_bin] = false
		options[:y2_bin] = false
	end
end

function ddsmMoveUp(sim::Simulation)
	@assert(sim.moveUpData.useMoveUp)
	
	# shorthand:
	ambulances = sim.ambulances
	stations = sim.stations
	numStations = sim.numStations
	currentTime = sim.time
	demand = sim.demand
	ddsmd = sim.moveUpData.ddsmData
	coverTimes = ddsmd.coverTimes
	coverTimeDemandPriorities = ddsmd.coverTimeDemandPriorities
	alpha = ddsmd.alpha
	slackWeight = ddsmd.slackWeight
	options = ddsmd.options
	
	# get movable ambulances (movableAmbs)
	ambMovable = [isAmbAvailableForMoveUp(amb) for amb in ambulances]
	movableAmbs = ambulances[ambMovable]
	numMovableAmbs = length(movableAmbs)
	
	if numMovableAmbs == 0
		return moveUpNull()
	end
	
	# calculate travel time for each available ambulance
	ambToStationTimes = zeros(Float, numMovableAmbs, numStations) # ambToStationTimes[i,j] = time for movable ambulance i to travel to station j
	for i = 1:numMovableAmbs
		ambToStationTimes[i,:] = ambMoveUpTravelTimes!(sim, movableAmbs[i])
	end
	ambToStationCost = ambToStationTimes * ddsmd.travelTimeCost
	
	# get demand point coverage data for the two target cover times
	# for coverTimes[ti], point j has demand pointDemands[ti][j] and is covered by stations pointStations[ti][j]
	pointStations = [] # pointStations[ti][j] is coverTimes[ti]
	pointDemands = []
	numPoints = Int[] # can have different number of points (really point sets) for each cover time
	demandPriorityArrivalRates = [getDemandMode!(demand, demandPriority, currentTime).arrivalRate for demandPriority in priorities]
	for demandPriority in coverTimeDemandPriorities
		# get demand point coverage data
		pointsCoverageMode = getPointsCoverageMode!(sim, demandPriority, currentTime)
		demandMode = getDemandMode!(demand, demandPriority, currentTime)
		pointSetsDemands = getPointSetsDemands!(sim, demandPriority, currentTime; pointsCoverageMode = pointsCoverageMode) * demandMode.rasterMultiplier
		pointSetsDemands *= sum(demandPriorityArrivalRates) / demandPriorityArrivalRates[Int(demandPriority)] # scale demand to be for all priorities, not just current demandPriority
		
		push!(pointStations, pointsCoverageMode.stationSets)
		push!(pointDemands, pointSetsDemands)
		push!(numPoints, length(pointSetsDemands))
	end
	# could have it so that pointStations, pointDemands, and numPoints are only updated if the demandSet has changed since last query,
	# but this may only shave off a couple of seconds per 100 sim days, which is not the current bottleneck
	
	######################
	# IP
	
	# shorthand
	a = numMovableAmbs
	np = numPoints # np[ti] is number of points for coverTimes[i]
	s = numStations
	t = coverTimes
	
	# using JuMP
	model = Model()
	jump_ge_0_19 = pkgVersions["JuMP"] >= v"0.19"
	solver = options[:solver]
	if jump_ge_0_19
		if solver == "cbc" set_optimizer(model, with_optimizer(Cbc.Optimizer, logLevel=0))
		elseif solver == "glpk" set_optimizer(model, with_optimizer(GLPK.Optimizer))
		elseif solver == "gurobi" @stdout_silent(set_optimizer(model, with_optimizer(Gurobi.Optimizer, OutputFlag=0)))
		end
	else
		if solver == "cbc" setsolver(model, CbcSolver())
		elseif solver == "glpk" setsolver(model, GLPKSolverMIP(presolve=true))
		elseif solver == "gurobi" @stdout_silent(setsolver(model, GurobiSolver(OutputFlag=0)))
		end
	end
	
	m = model # shorthand
	options[:x_bin] ? @variable(m, x[i=1:a,j=1:s], Bin) : @variable(m, 0 <= x[i=1:a,j=1:s] <= 1)
	options[:y11_bin] ? @variable(m, y11[p=1:np[1]], Bin) : @variable(m, 0 <= y11[p=1:np[1]] <= 1)
	options[:y12_bin] ? @variable(m, y12[p=1:np[1]], Bin) : @variable(m, 0 <= y12[p=1:np[1]] <= 1)
	options[:y2_bin] ? @variable(m, y2[p=1:np[2]], Bin) : @variable(m, 0 <= y2[p=1:np[2]] <= 1)
	# x[i,j] = 1 if ambulance: movableAmbs[i] should be moved to station: stations[j], 0 otherwise
	# y11[p,k] = 1 if demand point p is covered at least once within t[1], 0 otherwise
	# y12[p,k] = 1 if demand point p is covered at least twice within t[1], 0 otherwise
	# y2[p] = 1 if demand point p is covered at least once within t[2], 0 otherwise
	
	@variables(model, begin
		(s1 >= 0) # slack
		(s2 >= 0) # slack
	end)
	
	@constraints(model, begin
		(ambAtOneLocation[i=1:a], sum(x[i,:]) == 1) # each ambulance must be assigned to one station
		(pointCoverOrderY1[p=1:np[1]], y11[p] >= y12[p]) # single coverage before double coverage; not needed for y2
		(demandCoveredOnceT1, sum(y11[p] * pointDemands[1][p] for p=1:np[1]) + s1 >= alpha * sum(pointDemands[1])) # fraction (alpha) of demand covered once within t[1]
		(demandCoveredOnceT2, sum(y2[p] * pointDemands[2][p] for p=1:np[2]) + s2 >= sum(pointDemands[2])) # all demand covered once within t[2]
	end)
	
	if options[:z_var]
		@variable(model, z[j=1:s], Int) # z[j] = number of ambulances to assign to station j
		@constraints(model, begin
			(ambStationCount[j=1:s], z[j] == sum(x[:,j]))
			(pointCoverCountY1[p=1:np[1]], y11[p] + y12[p] <= sum(z[pointStations[1][p]]))
			(pointCoverCountY2[p=1:np[2]], sum(y2[p]) <= sum(z[pointStations[2][p]]))
		end)
	else
		@constraints(model, begin
			(pointCoverCountY1[p=1:np[1]], y11[p] + y12[p] <= sum(x[:,pointStations[1][p]]))
			(pointCoverCountY2[p=1:np[2]], sum(y2[p]) <= sum(x[:,pointStations[2][p]]))
		end)
	end
	
	@expressions(model, begin
		demandCoveredTwiceT1, sum(y12[p] * pointDemands[1][p] for p=1:np[1])
		totalAmbTravelCost, sum(x[i,j] * ambToStationCost[i,j] for i=1:a, j=1:s)
		slackCost, (s1 + s2) * slackWeight
	end)
	
	# solve
	if jump_ge_0_19
		@objective(model, Max, demandCoveredTwiceT1 - totalAmbTravelCost - slackCost)
		@stdout_silent optimize!(model)
		@assert(termination_status(model) == MOI.OPTIMAL)
	else
		@objective(model, :Max, demandCoveredTwiceT1 - totalAmbTravelCost - slackCost)
		status = @stdout_silent solve(model)
		@assert(status == :Optimal)
	end
	
	# get solution
	vals = Dict()
	vals[:x] = jump_ge_0_19 ? JuMP.value.(x) : getvalue(x)
	vals[:y11] = jump_ge_0_19 ? JuMP.value.(y11) : getvalue(y11)
	vals[:y12] = jump_ge_0_19 ? JuMP.value.(y12) : getvalue(y12)
	vals[:y2] = jump_ge_0_19 ? JuMP.value.(y2) : getvalue(y2)
	
	sol = convert(Array{Bool,2}, round.(vals[:x]))
	
	if checkMode
		@assert(all(sum(sol, dims=2) .== 1)) # check constraint: ambAtOneLocation
		# check that values are binary/integer
		for s in [:x, :y11, :y12, :y2]
			err = maximum(abs.(vals[s] - round.(vals[s])))
			@assert(err <= 10*eps(Float32), (s, err))
		end
	end
	
	ambStations = [stations[findfirst(sol[i,:])] for i=1:a]
	
	return movableAmbs, ambStations
end
