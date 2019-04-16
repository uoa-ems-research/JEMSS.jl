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
	coverFractionTargetT1::Float = 0.5, travelTimeCost::Float = 50.0, slackWeight::Float = 1000.0,
	coverTimeDemandPriorities::Vector{Priority} = [highPriority, lowPriority],
	options::Dict{Symbol,Any} = Dict{Symbol,Any}())
	# coverTimeDemandPriorities[i] is demand priority for coverTimes[i], where coverTimes[1] and [2] are the targets for ddsm
	
	@assert(0 <= coverFractionTargetT1 <= 1)
	@assert(length(coverTimeDemandPriorities) == 2)
	
	# initialise demand and demand coverage data if not already initialised
	sim.demand.initialised || initDemand!(sim)
	sim.demandCoverage.initialised || initDemandCoverage!(sim)
	demand = sim.demand # shorthand
	
	@assert(all([demandMode.rasterIndex for demandMode in demand.modes[demand.modeLookup[i,:]]] |> unique |> length == 1 for i = 1:demand.numSets)) # distribution of demand must be same for all demand priorities, in any given time period
	
	# two target cover times for ddsm
	coverTimes = [sim.demandCoverage.coverTimes[p] for p in coverTimeDemandPriorities]
	@assert(length(coverTimes) == 2)
	@assert(coverTimes[1] < coverTimes[2]) # change to warning?
	
	# set options
	ddsmd = sim.moveUpData.ddsmData # shorthand
	ddsmd.options[:solver] = "cbc" # can be slower than glpk, but more reliable for some reason
	ddsmd.options[:solver_args] = []
	ddsmd.options[:solver_kwargs] = []
	ddsmd.options[:v] = v"1"
	ddsmd.options[:z_var] = true
	merge!(ddsmd.options, Dict([:x_bin => true, :y11_bin => true, :y12_bin => true, :y2_bin => true]))
	ddsmOptions!(sim, options)
	
	ddsmd.coverFractionTargetT1 = coverFractionTargetT1
	ddsmd.travelTimeCost = travelTimeCost
	ddsmd.slackWeight = slackWeight
	ddsmd.coverTimeDemandPriorities = coverTimeDemandPriorities # coverTimeDemandPriorities[i] is demand priority for coverTimes[i]
	ddsmd.coverTimes = coverTimes
end

# change the options for ddsm
function ddsmOptions!(sim, options::Dict{Symbol,T}) where T <: Any
	options = merge!(sim.moveUpData.ddsmData.options, options) # update options
	
	@assert(in(options[:solver], ["cbc", "glpk", "gurobi"]))
	@assert(typeof(options[:v]) == VersionNumber)
	@assert(typeof(options[:z_var]) == Bool)
	
	if options[:solver] == "gurobi" try Gurobi; catch; options[:solver] = "cbc"; @warn("Failed to use Gurobi, using Cbc instead.") end end
	
	# options[:v] == v"0" # do nothing, values should already be set in options if using this
	options[:v] == v"1" && merge!(options, Dict([:x_bin => true, :y11_bin => true, :y12_bin => true, :y2_bin => true]))
	options[:v] == v"2" && merge!(options, Dict([:x_bin => true, :y11_bin => true, :y12_bin => false, :y2_bin => false]))
	options[:v] == v"3" && merge!(options, Dict([:x_bin => false, :y11_bin => true, :y12_bin => false, :y2_bin => false]))
	options[:v] == v"4" && merge!(options, Dict([:x_bin => false, :y11_bin => true, :y12_bin => true, :y2_bin => true]))
	if options[:y11_bin] == false && in(options[:solver], ["cbc", "glpk"]) @warn("Removing binary constraint for `y11` may not work with CBC or GLPK.") end
	if options[:x_bin] == false && options[:solver] == "gurobi" @warn("Removing binary constraint for `x` may not work with Gurobi.") end
	@assert(all(key -> haskey(options, key), [:x_bin, :y11_bin, :y12_bin, :y2_bin]))
	
	if isdefined(sim, :backup) sim.backup.moveUpData.ddsmData.options = options end # to keep options if sim is reset
	
	return options
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
	coverFractionTargetT1 = ddsmd.coverFractionTargetT1
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
	ambToStationCosts = ambToStationTimes * ddsmd.travelTimeCost
	
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
	@assert(isapprox(sum(pointDemands[1]), sum(pointDemands[2])))
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
	solver = options[:solver] # shorthand
	args = options[:solver_args] # shorthand
	kwargs = options[:solver_kwargs] # shorthand
	if jump_ge_0_19
		if solver == "cbc" set_optimizer(model, with_optimizer(Cbc.Optimizer, logLevel=0, args...; kwargs...))
		elseif solver == "glpk" set_optimizer(model, with_optimizer(GLPK.Optimizer, args...; kwargs...))
		elseif solver == "gurobi" @stdout_silent(set_optimizer(model, with_optimizer(Gurobi.Optimizer, OutputFlag=0, args...; kwargs...)))
		end
	else
		if solver == "cbc" setsolver(model, CbcSolver(args...; kwargs...))
		elseif solver == "glpk" setsolver(model, GLPKSolverMIP(args...; kwargs...))
		elseif solver == "gurobi" @stdout_silent(setsolver(model, GurobiSolver(OutputFlag=0, args...; kwargs...)))
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
		(demandCoveredOnceT1, sum(y11[p] * pointDemands[1][p] for p=1:np[1]) + s1 >= coverFractionTargetT1 * sum(pointDemands[1])) # fraction of demand covered once within t[1]
		(demandCoveredOnceT2, sum(y2[p] * pointDemands[2][p] for p=1:np[2]) + s2 >= sum(pointDemands[2])) # all demand covered once within t[2]
	end)
	
	if options[:z_var]
		@variable(model, z[j=1:s], Int) # z[j] = number of ambulances to assign to station j
		@constraints(model, begin
			(ambStationCount[j=1:s], z[j] == sum(x[:,j]))
			(pointCoverCountY1[p=1:np[1]], y11[p] + y12[p] <= sum(z[pointStations[1][p]]))
			(pointCoverCountY2[p=1:np[2]], y2[p] <= sum(z[pointStations[2][p]]))
		end)
	else
		@constraints(model, begin
			(pointCoverCountY1[p=1:np[1]], y11[p] + y12[p] <= sum(x[:,pointStations[1][p]]))
			(pointCoverCountY2[p=1:np[2]], y2[p] <= sum(x[:,pointStations[2][p]]))
		end)
	end
	
	@expressions(model, begin
		demandCoveredTwiceT1, sum(y12[p] * pointDemands[1][p] for p=1:np[1])
		totalAmbTravelCost, sum(x[i,j] * ambToStationCosts[i,j] for i=1:a, j=1:s)
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
		for sym in [:x, :y11, :y12, :y2]
			err = maximum(abs.(vals[sym] - round.(vals[sym])))
			@assert(err <= 10*eps(Float32), (sym, err))
		end
	end
	
	ambStations = [stations[findfirst(sol[i,:])] for i=1:a]
	
	return movableAmbs, ambStations
end
