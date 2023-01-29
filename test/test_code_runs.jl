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

# tests to check that code runs without crashing, but not checking that output is expected

# check that different sim configs can be opened, and sim can run
@testset "sim configs" begin
    @assert(isdir("data/cities/small/1/generated"))
    simConfigFolder = "data/cities/small/1/sim_configs"
    for configFilename in readdir(simConfigFolder)
        filename = joinpath(pwd(), simConfigFolder, configFilename)
        sim = initSim(filename)
        simulate!(sim)
        @test true
    end
end

# test the scripts in /example/scripts
if true
    @info("Skipped testing of example scripts") # default is to not run this test set
else
    @testset "example scripts" begin
        scriptsFolder = joinpath(dirname(pathof(JEMSS)), "..", "example", "scripts")

        # local search scripts, tested with a very small sim for speed
        cd("data/cities/small/6") do
            runGenConfig("gen_config.xml", overwriteOutputPath=true, doPrint=false)
            isdir("output") || mkdir("output")
            for script in ("deployment_local_search.jl", "deployment_search.jl", "deployment_test.jl",
                "priority_list_local_search.jl", "nested_comp_table_local_search.jl")
                @info(string("Running script: ", script))
                include(joinpath(scriptsFolder, script))
                @test true
                println()
            end
            rm("output", recursive=true)
        end

        # transient.jl
        cd("data/cities/small/2") do
            runGenConfig("gen_config.xml", overwriteOutputPath=true, doPrint=false)
            isdir("output") || mkdir("output")
            @info(string("Running script: transient.jl"))
            include(joinpath(scriptsFolder, "transient.jl"))
            @test true
            println()
            rm("output", recursive=true)
        end
    end
end

@testset "generate calls" begin
    cd("data/cities/small/3") do
        runGenConfig("gen_config_calls.xml", overwriteOutputPath=true, doPrint=false) # single calls file, limit calls by count
        runGenConfig("gen_config_calls_2.xml", overwriteOutputPath=true, doPrint=false) # single calls file, with maximum arrival time for last call
        runGenConfig("gen_config_calls_3.xml", overwriteOutputPath=true, doPrint=false) # single calls file, with minimum arrival time for last call
        genConfig = runGenConfig("gen_config_calls_multiple.xml", overwriteOutputPath=true, doPrint=false) # multiple calls files
        callSets = runGenConfigCalls(genConfig, doPrint=false, writeFile=false)
        @assert(isa(callSets, Vector{Vector{Call}}) && length(callSets) == genConfig.numCallsFiles)
        @test true
    end
end

@testset "set sim calls" begin
    sim = initSim("data/cities/small/1/sim_config.xml")
    genConfig = readGenConfig("data/cities/small/1/gen_config.xml")
    callSets = [makeCalls(genConfig) for i = 1:3]
    for calls in callSets
        for call in calls
            call.index = 0
        end # test that setSimCalls! can handle this
        @assert(calls[1].nearestNodeIndex == nullIndex) # test that setSimCalls! can handle this
        setSimCalls!(sim, calls)
        @assert(sim.calls[1].index == 1)
        simulate!(sim)
    end
    @test true
end

@testset "sim replications" begin
    cd("data/cities/small/5") do
        runGenConfig("gen_config.xml", overwriteOutputPath=true, doPrint=false)
        sim = initSim("sim_config.xml")
        @assert(sim.numReps == 3)
        simulateReps!(sim.reps)
        @assert(all(rep -> rep.complete, sim.reps))
        writeStatsFiles(sim, sim.reps)
        resetReps!(sim.reps)
        @test true
    end
end

@testset "sim replications parallel" begin
    cd("data/cities/small/5") do
        runGenConfig("gen_config.xml", overwriteOutputPath=true, doPrint=false)
        sim = initSim("sim_config.xml")
        @assert(sim.numReps == 3)
        simulateReps!(sim.reps; parallel=true, numThreads=2)
        @assert(all(rep -> rep.complete, sim.reps))
        writeStatsFiles(sim, sim.reps)
        resetReps!(sim.reps; parallel=true, numThreads=2)
        @test true
    end
end

# test that mexclp can be solved, but not solution correctness
@testset "mexclp" begin
    @assert(isdir("data/cities/small/1/generated"))
    filename = joinpath(pwd(), "data/cities/small/1/mexclp/sim_config.xml")
    sim = initSim(filename)
    # solve mexclp with various kwarg values
    demandWeights = Dict([p => 1.0 for p in priorities])
    demandWeights[lowPriority] = 0
    stationsNumAmbs, deployment = solveMexclp!(sim; busyFraction=0.4, demandWeights=demandWeights, stationCapacities=[2 for i = 1:sim.numStations])
    solveMexclp!(sim; busyFraction=-0.4) # test formulation where adding an ambulance at a station can reduce coverage
    # check that solution can be applied
    applyStationsNumAmbs!(sim, stationsNumAmbs)
    applyDeployment!(sim, deployment)
    @test true
end

@testset "cover bound" begin
    runGenConfig("data/cities/small/7/gen_config.xml", overwriteOutputPath=true, doPrint=false)
    sim = initSim("data/cities/small/7/sim_config.xml")
    genConfig = readGenConfig(sim.inputFiles["callGenConfig"].path)
    coverBoundSim = initCoverBoundSim(numReps=30, numAmbs=sim.numAmbs, warmUpDuration=0.1, minLastCallArrivalTime=0.1 + 10,
        interarrivalTimeDistrRng=genConfig.interarrivalTimeDistrRng, ambBusyDurationSeed=1)
    ambBusyDurationsToSample = [1:60;] / 60 / 24 # 1 minute intervals, up to 1 hour
    queuedDurationsToSample = [0:maximum(sim.targetResponseDurations)*24*60;] / 60 / 24 # 1 minute intervals
    coverBound = calcCoverBound!(sim; coverBoundSim=coverBoundSim, ambBusyDurationsToSample=ambBusyDurationsToSample, queuedDurationsToSample=queuedDurationsToSample)
    calcCoverBound!(coverBound; accountForQueuedDurations=false) # also calculate bound without accounting for call queueing durations (original cover bound)
    # @show coverBound.sim.bound # cover bound result
    simulateCoverBoundLowerBound!(coverBound) # also simulate a lower bound on the cover bound
    @test true
end
