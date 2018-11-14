# deployment (allocation) of ambulances to stations

# make a single random deployment
# ignores station capacities, assumes ambulances are homogeneous
function makeRandDeployment(numAmbs::Int, numStations::Int;
	rng::AbstractRNG = Base.GLOBAL_RNG)::Deployment
	return sort(rand(rng, 1:numStations, numAmbs))
end

# generate a number of unique random deployments
# ignores station capacities, assumes ambulances are homogeneous
function makeRandDeployments(numAmbs::Int, numStations::Int, numDeployments::Int;
	rng::AbstractRNG = Base.GLOBAL_RNG)
	@assert(numAmbs >= 1)
	@assert(numStations >= 1)
	@assert(numDeployments >= 1)
	@assert(numDeployments <= binomial(numAmbs + numStations - 1, numStations - 1))
	
	deployments = Vector{Deployment}() # deployments, deployments[i][j] gives station index that ambulance j should be deployed to, for ith deployment
	while length(deployments) < numDeployments
		newDeployment = makeRandDeployment(numAmbs, numStations; rng = rng)
		if !in(newDeployment, deployments)
			push!(deployments, newDeployment)
		end
	end
	
	return deployments
end
# function makeRandDeployments(sim::Simulation, numDeployments::Int)
	# return makeRandDeployments(sim.numAmbs, sim.numStations, numDeployments)
# end

# set ambulance to be located at given station
# this is hacky, expects ambulance will only need changing, this may not always be true
# mutates: ambulance
function setAmbStation!(ambulance::Ambulance, station::Station)
	ambulance.stationIndex = station.index
	ambulance.route.endLoc = station.location
	ambulance.route.endFNode = station.nearestNodeIndex
end

# apply deployment 'deployment' for sim.ambulances and sim.stations
# mutates: sim.ambulances
function applyDeployment!(sim::Simulation, deployment::Deployment)
	n = length(deployment)
	@assert(n == sim.numAmbs)
	@assert(all(d -> in(d, 1:sim.numStations), deployment))
	# @assert(all(1 .<= deployment .<= sim.numStations))
	for i = 1:n
		setAmbStation!(sim.ambulances[i], sim.stations[deployment[i]])
	end
end

# runs simulation for the deployment
# mutates: sim
function simulateDeployment!(sim::Simulation, deployment::Deployment)
	resetSim!(sim)
	applyDeployment!(sim, deployment)
	simulateToEnd!(sim)
end

# runs simulation for each deployment
# returns vector of function applied to end of each simulation run
# mutates: sim
function simulateDeployments!(sim::Simulation, deployments::Vector{Deployment}, f::Function;
	showEta::Bool = false)
	numDeployments = length(deployments)
	results = Vector{Any}(numDeployments)
	t = time() # for estimating expected time remaining
	for i = 1:numDeployments
		simulateDeployment!(sim, deployments[i])
		results[i] = f(sim)
		if showEta
			remTime = (time() - t) / i * (numDeployments - i) / 60 # minutes
			print(rpad(string(" remaining run time (minutes): ", ceil(remTime,1)), 100, " "), "\r")
		end
	end
	resetSim!(sim)
	return results
end
