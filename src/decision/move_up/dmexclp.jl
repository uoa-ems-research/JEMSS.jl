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

# Dynamic Maximum Expected Coverage Location Problem (DMEXCLP)
# From paper: "An efficient heuristic for real-time ambulance redeployment"

"""
	function initDmexclp!(sim::Simulation;
		busyFraction::Float = 0.5,
		demandWeights::Dict{Priority,Float} = Dict([p => 1.0 for p in priorities]))
Initialise data for the Dynamic Maximum Expected Coverage Location Problem (DMEXCLP).

# Keyword arguments
- `busyFraction` is the fraction of time that ambulances are busy; should be within [0,1] though this is not enforced
- `demandWeights` is the weight of each demand priority on the objective function

Mutates: `sim.moveUpData.dmexclpData`
"""
function initDmexclp!(sim::Simulation;
	busyFraction::Float = 0.5,
	demandWeights::Dict{Priority,Float} = Dict([p => 1.0 for p in priorities]))
	
	# initialise demand and demand coverage data if not already initialised
	sim.demand.initialised || initDemand!(sim)
	sim.demandCoverage.initialised || initDemandCoverage!(sim)
	
	# shorthand
	numAmbs = sim.numAmbs
	numStations = sim.numStations
	
	dcd = sim.moveUpData.dmexclpData # shorthand
	dcd.busyFraction = busyFraction
	dcd.demandWeights = demandWeights
	
	# calculate cover benefit values, for single demand
	dcd.marginalBenefit = (busyFraction.^[0:numAmbs-1;])*(1-busyFraction)
	
	# values that will be calculated when needed
	dcd.stationNumFreeAmbs = zeros(Int, numStations)
	dcd.stationMarginalCoverages = zeros(Float, numStations) # stationMarginalCoverages[i] gives extra coverage provided from placing newly freed ambulance at station i
	# dcd.pointSetsCoverCounts = [zeros(Int, length(pointsCoverageMode.pointSets)) for pointsCoverageMode in sim.demandCoverage.pointsCoverageModes] # pointSetsCoverCounts[i][j] = number of free ambulances covering node set j, for demand.pointsCoverageModes i
end

function dmexclpMoveUp(sim::Simulation, newlyFreedAmb::Ambulance)
	@assert(sim.moveUpData.useMoveUp)
	@assert(newlyFreedAmb.status != ambGoingToCall)
	@assert(sim.demand.initialised && sim.demandCoverage.initialised)
	
	# shorthand names:
	dcd = sim.moveUpData.dmexclpData
	ambulances = sim.ambulances
	stations = sim.stations
	numStations = sim.numStations
	
	# calculate the number of free ambulances at (or travelling to) each station
	dcd.stationNumFreeAmbs[:] .= 0
	for (i,ambulance) in enumerate(ambulances)
		# do not count newly freed ambulance, it has not been assigned a station
		if isAmbAvailableForMoveUp(ambulance) && i != newlyFreedAmb.index
			dcd.stationNumFreeAmbs[ambulance.stationIndex] += 1
		end
	end
	
	# ignoring newly freed amb, count number of ambulances covering each point set
	demandsPointSetsCoverCounts = calcPointSetsCoverCounts!(sim, sim.time, dcd.stationNumFreeAmbs)
	
	# find station allocation for newly freed ambulance that gives greatest
	# increase in expected demand coverage
	dcd.stationMarginalCoverages[:] .= 0.0
	for demandPriority in priorities
		if !haskey(dcd.demandWeights, demandPriority) || dcd.demandWeights[demandPriority] == 0
			continue
		end
		pointSetsCoverCounts = demandsPointSetsCoverCounts[demandPriority]
		pointsCoverageMode = getPointsCoverageMode!(sim, demandPriority, sim.time)
		demandMode = getDemandMode!(sim.demand, demandPriority, sim.time)
		pointSetsDemands = sim.demandCoverage.pointSetsDemands[pointsCoverageMode.index, demandMode.rasterIndex]
		# to do: if marginal coverage has already been calculated for the combination of pointsCoverageMode and rasterIndex (but for a different demand priority), the marginal coverage value should be reused (accounting for differences in old and new values of demandMode.rasterMultiplier and dcd.demandWeights[demandPriority]) to save on computation.
		for i = 1:length(pointsCoverageMode.pointSets)
			pointSetDemand = pointSetsDemands[i] * demandMode.rasterMultiplier
			pointSetMarginalCoverage = pointSetDemand * dcd.marginalBenefit[pointSetsCoverCounts[i]+1] * dcd.demandWeights[demandPriority]
			for j in pointsCoverageMode.stationSets[i]
				dcd.stationMarginalCoverages[j] += pointSetMarginalCoverage
			end
		end
	end
	(bestMarginalCoverage, bestStationIndex) = findmax(dcd.stationMarginalCoverages)
	
	return [newlyFreedAmb], [stations[bestStationIndex]]
end
