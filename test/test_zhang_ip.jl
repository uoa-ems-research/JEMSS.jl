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

@testset "zhang ip" begin
	numStations = 2
	stationSlots = [1, 1, 2, 2] # each station can hold two ambulances
	benefitSlots = [4.0, 1.0, 3.0, 2.0] # benefitSlots[i] gives value of placing ambulance at stationSlots[i]
	@assert(all(j -> issorted(benefitSlots[stationSlots .== j], lt=<=, rev=true), numStations)) # value of placing kth amb at station j > value of k+1th amb
	
	function solveZhangIpTest(stationSlots, benefitSlots, ambToStationCosts, numMovableAmbs, numAtHospitalAmbs)
		if typeof(ambToStationCosts) <: Vector
			ambToStationCosts = convert(Array{Float,2}, ambToStationCosts')
		end
		solution1 = JEMSS.solveZhangIp(stationSlots, benefitSlots, ambToStationCosts, numMovableAmbs, numAtHospitalAmbs, true)
		solution2 = JEMSS.solveZhangIpAssignmentProblem(stationSlots, benefitSlots, ambToStationCosts, numMovableAmbs, numAtHospitalAmbs)
		return solution1 == solution2 ? solution1 : (solution1, solution2)
	end
	
	# case: no ambs
	ambToStationCosts = zeros(Float, 0, numStations) # ambToStationCosts[i,j] is cost of assigning amb i to station j (excludes benefit)
	expectedSolution = (Int[], Int[])
	@test expectedSolution == solveZhangIpTest(stationSlots, benefitSlots, ambToStationCosts, 0, 0)
	
	## 1 amb cases
	
	# case: 1 amb at station (not moved)
	ambToStationCosts = Float[1.5 0] # amb at station 2; cost to move is greater than benefit increase
	expectedSolution = ([2], Int[])
	@test expectedSolution == solveZhangIpTest(stationSlots, benefitSlots, ambToStationCosts, 1, 0)
	
	# case: 1 amb at station (moved)
	ambToStationCosts = Float[0.5 0] # amb at station 2; cost to move is less than benefit increase
	expectedSolution = ([1], Int[]) # moved to station 1
	@test expectedSolution == solveZhangIpTest(stationSlots, benefitSlots, ambToStationCosts, 1, 0)
	
	# case: 1 amb at station (assigned to station even with negative objective value)
	ambToStationCosts = Float[5 5]
	expectedSolution = ([1], Int[])
	@test expectedSolution == solveZhangIpTest(stationSlots, benefitSlots, ambToStationCosts, 1, 0)
	
	# case: one amb at hospital (assigned to station 1)
	ambToStationCosts = Float[1.5 1] # value for station 1: 4 - 1.5 = 3.5; for station 2: 3 - 1 = 2
	expectedSolution = (Int[], [1])
	@test expectedSolution == solveZhangIpTest(stationSlots, benefitSlots, ambToStationCosts, 0, 1)
	
	# case: one amb at hospital (not assigned to a station)
	ambToStationCosts = Float[5 5] # costs outweight benefits, will not redeploy amb
	expectedSolution = (Int[], [nullIndex])
	@test expectedSolution == solveZhangIpTest(stationSlots, benefitSlots, ambToStationCosts, 0, 1)
	
	## 2 amb cases
	
	# case: 2 ambs at stations (not moved)
	ambToStationCosts = Float[0 1; 1 0]
	expectedSolution = ([1,2], Int[])
	@test expectedSolution == solveZhangIpTest(stationSlots, benefitSlots, ambToStationCosts, 2, 0)
	
	# case: 2 ambs at 1 station (not moved)
	ambToStationCosts = Float[0 4; 0 4]
	expectedSolution = ([1,1], Int[])
	@test expectedSolution == solveZhangIpTest(stationSlots, benefitSlots, ambToStationCosts, 2, 0)
	
	# case: 2 ambs at 1 station (amb 2 moved)
	ambToStationCosts = Float[0 2; 0 1]
	expectedSolution = ([1,2], Int[])
	@test expectedSolution == solveZhangIpTest(stationSlots, benefitSlots, ambToStationCosts, 2, 0)
	
	# case: 1 amb at station, 1 at hospital; amb at hospital keeps amb from moving to better station
	ambToStationCosts = Float[0.5 0; 1 1] # amb 1 would move to station 1 if second amb wasn't at hospital
	# note: ambToStationCosts[1,1] + ambToStationCosts[2,2] > ambToStationCosts[1,2] + ambToStationCosts[2,1]
	@test ([1], Int[]) == solveZhangIpTest(stationSlots, benefitSlots, ambToStationCosts[1,:], 1, 0)
	@test ([2], [1]) == solveZhangIpTest(stationSlots, benefitSlots, ambToStationCosts, 1, 1)
	
	# case: 1 amb at station, 1 at hospital; amb at hospital causes amb to move to worse station
	ambToStationCosts = Float[0 0.5; 1 2]
	# note: ambToStationCosts[1,1] + ambToStationCosts[2,2] > ambToStationCosts[1,2] + ambToStationCosts[2,1]
	@test ([1], Int[]) == solveZhangIpTest(stationSlots, benefitSlots, ambToStationCosts[1,:], 1, 0) # without amb at hospital
	@test ([2], [1]) == solveZhangIpTest(stationSlots, benefitSlots, ambToStationCosts, 1, 1)
	
	## 4 amb cases
	
	# case: 4 ambs at stations
	ambToStationCosts = Float[0 1; 0 1; 1 0; 1 0]
	@test ([1,1,2,2], Int[]) == solveZhangIpTest(stationSlots, benefitSlots, ambToStationCosts, 4, 0)
	
	# case: 3 ambs at stations, 1 at hospital; amb at hospital causes amb 3 to move from station 2 to 1
	ambToStationCosts = Float[0 2; 1 0; 0.5 0; 1 0.1]
	# note: ambToStationCosts[3,1] + ambToStationCosts[4,2] < ambToStationCosts[3,2] + ambToStationCosts[4,1]
	@test ([1,2,2], Int[]) == solveZhangIpTest(stationSlots, benefitSlots, ambToStationCosts[1:3,:], 3, 0) # without amb at hospital
	@test ([1,2,1], [2]) == solveZhangIpTest(stationSlots, benefitSlots, ambToStationCosts, 3, 1) # with amb at hospital
	
	# case: 2 ambs at station 1, 2 at hospital; ambs at hospital cause ambs at station 1 to move to station 2
	ambToStationCosts = Float[0 2.2; 0 2.1; -2 4; -2 4] # ambs at station 1 would not move if not for ambs at hospitals
	@test ([1,1], Int[]) == solveZhangIpTest(stationSlots, benefitSlots, ambToStationCosts[1:2,:], 2, 0) # no ambs at hospital
	@test ([1,2], [1]) == solveZhangIpTest(stationSlots, benefitSlots, ambToStationCosts[1:3,:], 2, 1) # 1 amb at hospital
	@test ([2,2], [1,1]) == solveZhangIpTest(stationSlots, benefitSlots, ambToStationCosts, 2, 2) # 2 ambs at hospital
	
end
