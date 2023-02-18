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

@testset "get travel mode" begin
    sim = initSim("data/cities/small/1/sim_config.xml") # see travel.csv file in same folder
    travel = sim.travel # shorthand

    # first test getTravelSetsStartTimesIndex(), which getTravelMode!() relies on
    f(t::Float) = JEMSS.getTravelSetsStartTimesIndex(travel, t)
    result = true
    travelSetsStartTimes = [0, 0.2, 0.5, 0.8]
    for t = 0:0.01:1.0
        expectedValue = findlast(t .>= travelSetsStartTimes)
        result &= f(t) == expectedValue
    end
    @test result

    # test getTravelMode!()
    result = true
    f(p::Priority, t::Float) = JEMSS.getTravelMode!(travel, p, sim.time; startTime=t).index
    expectedValues = Dict(
        (highPriority, 0.0) => 1,
        (medPriority, 0.0) => 2,
        (lowPriority, 0.0) => 3,
        (highPriority, 1.0) => 2,
        (medPriority, 1.0) => 3,
        (lowPriority, 1.0) => 4,
    )
    for priority in priorities, t = (0.0, 1.0)
        result &= f(priority, t) == expectedValues[(priority, t)]
    end
    @test result

    # check that getTravelMode!() updates travel.recentSetsStartTimesIndex
    function g(t::Float)
        JEMSS.getTravelMode!(travel, highPriority, t)
        return travel.recentSetsStartTimesIndex
    end
    result = true
    travelSetsStartTimes = [0, 0.2, 0.5, 0.8]
    for t = 0:0.01:1.0
        expectedValue = findlast(t .>= travelSetsStartTimes)
        result &= g(t) == expectedValue
    end
    @test result
end