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

# For determining transient period duration at start and end of simulation, based on call response durations.
# Useful for deciding how duration of simulation warm-up and cool-down periods.
# Loads sim from config, splits calls into batches, runs simulation for each batch,
# then saves and plots the moving average of response durations.

using JEMSS
using Statistics

# parameters
configFilename = "sim_config.xml"
outputPath = "output"
displayPlots = false

println("Loading sim")
sim = initSim(configFilename, doPrint = false)

isdir(outputPath) || mkdir(outputPath)

callArrivalTimes = [c.arrivalTime for c in sim.calls] # backup, as call arrival times will change when calls are batched

# break calls into call batches, each with callBatchNumCalls many calls, or numCallBatches many call batches
# mutates: calls
function createCallBatchesByCount!(calls::Vector{Call}, startTime::Float;
	callBatchNumCalls::Int = nullIndex, numCallBatches::Int = nullIndex)
	
	@assert(all(c -> c.status == callNullStatus, calls))
	
	numCalls = length(calls)
	@assert(sum([callBatchNumCalls, numCallBatches] .!= nullIndex) == 1)
	if callBatchNumCalls != nullIndex
		numCallBatches = Int(floor(numCalls / callBatchNumCalls))
	elseif numCallBatches != nullIndex
		callBatchNumCalls = Int(floor(numCalls / numCallBatches))
	end
	
	callBatches = Vector{Vector{Call}}(undef,numCallBatches)
	callBatchStartTimes = Vector{Float}(undef,numCallBatches)
	callBatchFirstCallIndex = callBatchLastCallIndex = 0 # first and last call indices for current call batch
	for callBatchIndex = 1:numCallBatches
		# get call batch start time
		if callBatchIndex == 1
			callBatchStartTimes[1] = startTime
		else
			callBatchStartTimes[callBatchIndex] = calls[callBatchLastCallIndex].arrivalTime
		end
		# find first and last calls for call batch callBatchIndex
		callBatchFirstCallIndex = callBatchLastCallIndex + 1
		callBatchLastCallIndex = callBatchFirstCallIndex + callBatchNumCalls - 1
		callBatches[callBatchIndex] = calls[callBatchFirstCallIndex:callBatchLastCallIndex]
	end
	
	# renumber calls from 1 to n within each call batch, set each call batch to start at startTime
	for (callBatchIndex, callBatch) in enumerate(callBatches)
		for (callIndex, call) in enumerate(callBatch)
			call.index = callIndex
			call.arrivalTime -= (callBatchStartTimes[callBatchIndex] - startTime)
		end
	end
	
	return callBatches
end

# run simulation for each callBatch in callBatches
# mutates: sim, callBatches
function simulateCallBatches!(sim::Simulation, callBatches::Vector{Vector{Call}})
	# run simulation for each call batch
	numCallBatches = length(callBatches)
	for (callBatchIndex, callBatch) in enumerate(callBatches)
		print("\rSimulating call batch ", callBatchIndex, " of ", numCallBatches)
		setSimCalls!(sim, callBatch)
		simulate!(sim)
		sim.calls = Call[] # remove call batch
	end
	println()
end

function calcAvgResponseDurations(callBatches::Vector{Vector{Call}})
	# for callBatches that have already been used in simulation, calculate avgResponseDurations,
	# where avgResponseDurations[i] is the average response duration of call i from all batches
	@assert(all(callBatch -> all(call -> call.responseDuration != nullTime, callBatch), callBatches))
	avgResponseDurations = [mean([callBatches[i][j].responseDuration for i = 1:length(callBatches)]) for j = 1:length(callBatches[1])] * (24*60)
	return avgResponseDurations
end

function calcMovingAvgResponseDurations(avgResponseDurations::Vector{Float}, movingWindowSize::Int)
	callBatchNumCalls = length(avgResponseDurations)
	
	# # for constant window size:
	# numWindows = max(0, callBatchNumCalls - movingWindowSize + 1)
	# movingAvgResponseDurations = [mean(avgResponseDurations[i + (0:movingWindowSize-1)]) for i = 1:numWindows]
	
	# # for window size = min(window size, number of calls so far):
	# numWindows = callBatchNumCalls
	# movingAvgResponseDurations = [mean(avgResponseDurations[max(1,i-(movingWindowSize-1)):i]) for i = 1:numWindows]
	
	# for window size = 1 + 2 * min(movingWindowSize, max allowable window size):
	# (each average value is at centre of window)
	numWindows = callBatchNumCalls
	movingAvgResponseDurations = Vector{Float}(undef,numWindows)
	for windowIndex = 1:numWindows
		currentWindowSize = min(movingWindowSize, windowIndex - 1, numWindows - windowIndex)
		movingAvgResponseDurations[windowIndex] = mean(avgResponseDurations[windowIndex .+ (-currentWindowSize : currentWindowSize)])
	end
	
	# # for window size = 1 + 2 * min(movingWindowSize, max allowable window size):
	# # (each average value is at centre of window, last n (=movingWindowSize) average values are not calculated)
	# # this follows Welch's procedure in "On the problem of the initial transient in steady-state simulation"
	# numWindows = max(0, callBatchNumCalls - movingWindowSize)
	# movingAvgResponseDurations = Vector{Float}(undef,numWindows)
	# for windowIndex = 1:numWindows
		# currentWindowSize = min(movingWindowSize, max(0, windowIndex-1))
		# movingAvgResponseDurations[windowIndex] = mean(avgResponseDurations[windowIndex .+ (-currentWindowSize : currentWindowSize)])
	# end
	
	return movingAvgResponseDurations
end

# For calls in sim.calls, split into n batches (n = callSetsNumCallBatches[i], i=1,...), then simulate.
# The response duration of the jth call is averaged over all batches.
function getCallSetsAvgResponseDurations(sim, callSetsNumCallBatches::Vector{Int})::Vector{Vector{Float}}
	numCallSets = length(callSetsNumCallBatches)
	callSetsAvgResponseDurations = []
	calls = sim.calls # shorthand
	for (callSetIndex, numCallBatches) in enumerate(callSetsNumCallBatches)
		println("Call set ", callSetIndex, " of ", numCallSets)
		
		# create batches (for current set), simulate for each, calculate average response durations
		reset!(sim)
		callBatches = createCallBatchesByCount!(calls, sim.startTime; numCallBatches = numCallBatches)
		simulateCallBatches!(sim, callBatches)
		push!(callSetsAvgResponseDurations, calcAvgResponseDurations(callBatches))
		
		# undo changes from createCallBatchesByCount!(), and reset calls
		for (i,c) in enumerate(calls)
			c.index = i
			c.arrivalTime = callArrivalTimes[i]
		end
		setSimCalls!(sim, calls)
	end
	
	return callSetsAvgResponseDurations
end

function writeCallSetsAvgResponseDurationsFile(filename::String, callSetsAvgResponseDurations::Vector{Vector{Float}}, callSetsNumCallBatches::Vector{Int}, movingWindowSizes::Vector{Int})
	tables = Table[]
	for (i, numCallBatches) in enumerate(callSetsNumCallBatches)
		avgResponseDurations = callSetsAvgResponseDurations[i]
		m = length(avgResponseDurations)
		header = ["callIndex"]
		data = Array{Any,2}(undef, m, length(movingWindowSizes) + 1) # to be filled
		data[:] .= "" # empty values
		data[:,1] = 1:m
		for (j, movingWindowSize) in enumerate(movingWindowSizes)
			push!(header, "windowSize_$movingWindowSize")
			movingAvgResponseDurations = calcMovingAvgResponseDurations(avgResponseDurations, movingWindowSize)
			if movingAvgResponseDurations != []
				data[1:length(movingAvgResponseDurations), j+1] = movingAvgResponseDurations
			end
		end
		push!(tables, Table("numCallBatches_$numCallBatches", header, data))
	end
	writeTablesToFile(filename, tables)
end

# Given the response duration of the jth call averaged over all batches, moving averages are calculated using
# each value in movingWindowSizes. One plot is produced for each value in callSetsNumCallBatches.
function plotMovingAvgResponseDurations(callSetsAvgResponseDurations::Vector{Vector{Float}}, callSetsNumCallBatches::Vector{Int}, movingWindowSizes::Vector{Int})
	numCallSets = length(callSetsNumCallBatches)
	@assert(numCallSets == length(callSetsAvgResponseDurations))
	plots = []
	for i = 1:numCallSets
		# create plot of moving average response durations with different moving window sizes
		push!(plots, Plots.plot(title = string("Number of call batches: ", callSetsNumCallBatches[i]),
			xaxis = "Call index", yaxis = "Mean response duration (minutes)"))
		for movingWindowSize in movingWindowSizes
			movingAvgResponseDurations = calcMovingAvgResponseDurations(callSetsAvgResponseDurations[i], movingWindowSize)
			if movingAvgResponseDurations != []
				Plots.plot!(plots[i], movingAvgResponseDurations, label = string(movingWindowSize));
			end
		end
	end
	return plots
end

println("Simulating with different call batches")
callSetsNumCallBatches = [5,10] # for runtests.jl
# callSetsNumCallBatches = [100:100:1000;] # for actual use
movingWindowSizes = [0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000];
callSetsAvgResponseDurations = getCallSetsAvgResponseDurations(sim, callSetsNumCallBatches)

println("Saving results to file")
filename = joinpath(outputPath, "transient.csv")
writeCallSetsAvgResponseDurationsFile(filename, callSetsAvgResponseDurations, callSetsNumCallBatches, movingWindowSizes)

if displayPlots
	println("Plotting")
	import Plots
	Plots.plotly()
	plots = plotMovingAvgResponseDurations(callSetsAvgResponseDurations, callSetsNumCallBatches, movingWindowSizes)
	for plot in plots display(plot) end
end
