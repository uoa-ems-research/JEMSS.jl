# deployment (allocation) of ambulances to stations

# make a single random deployment policy
# ignores station capacities, assumes ambulances are homogeneous
function makeRandDeploymentPolicy(numAmbs::Int, numStations::Int)::Depol
	return sort(rand(1:numStations, numAmbs))
end

# generate a number of unique random deployment policies
# ignores station capacities, assumes ambulances are homogeneous
function makeRandDeploymentPolicies(numAmbs::Int, numStations::Int, numDepols::Int)
	assert(numAmbs >= 1)
	assert(numStations >= 1)
	assert(numDepols >= 1)
	assert(numDepols <= binomial(numAmbs + numStations - 1, numStations - 1))
	
	depols = Vector{Depol}() # deployment policies, depols[i][j] gives station index that ambulance j should be deployed to, for ith deployment policy
	while length(depols) < numDepols
		newDepol = makeRandDeploymentPolicy(numAmbs, numStations)
		if !in(newDepol, depols)
			push!(depols, newDepol)
		end
	end
	
	return depols
end
# function makeRandDeploymentPolicies(sim::Simulation, numDepols::Int)
	# return makeRandDeploymentPolicies(length(sim.ambulances), length(sim.stations), numDepols)
# end

# set ambulance to be located at given station
# this is hacky, expects ambulance will only need changing, this may not always be true
# mutates: ambulance
function setAmbStation!(ambulance::Ambulance, station::Station)
	ambulance.stationIndex = station.index
	ambulance.route.endLoc = station.location
	ambulance.route.endFNode = station.nearestNodeIndex
end

# apply deployment policy 'depol' for sim.ambulances and sim.stations
# mutates: sim.ambulances
function applyDeploymentPolicy!(sim::Simulation, depol::Depol)
	n = length(depol)
	assert(n == length(sim.ambulances))
	assert(all(d -> in(d, 1:length(sim.stations)), depol))
	# assert(all(1 .<= depol .<= length(sim.stations)))
	for i = 1:n
		setAmbStation!(sim.ambulances[i], sim.stations[depol[i]])
	end
end

# runs simulation for the deployment policy
# mutates: sim
function simulateDeploymentPolicy!(sim::Simulation, depol::Depol)
	resetSim!(sim)
	applyDeploymentPolicy!(sim, depol)
	simulateToEnd!(sim)
end

# runs simulation for each deployment policy
# returns vector of function applied to end of each simulation run
# mutates: sim
function simulateDeploymentPolicies!(sim::Simulation, depols::Vector{Depol}, f::Function;
	showEta::Bool = false)
	numDepols = length(depols)
	results = Vector{Any}(numDepols)
	t = time() # for estimating expected time remaining
	for i = 1:numDepols
		simulateDeploymentPolicy!(sim, depols[i])
		results[i] = f(sim)
		if showEta
			remTime = (time() - t) / i * (numDepols - i) / 60 # minutes
			print(rpad(string(" remaining run time (minutes): ", ceil(remTime,1)), 100, " "), "\r")
		end
	end
	resetSim!(sim)
	return results
end
