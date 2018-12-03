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


function getCallResponseTimes(sim::Simulation)
	@assert(sim.complete)
	@assert(all(call -> call.responseTime >= 0, sim.calls))
	return map(call -> call.responseTime, sim.calls)
end

function getAvgCallResponseTime(sim::Simulation; useMinutes::Bool = false)
	@assert(sim.complete)
	@assert(all(call -> call.responseTime >= 0, sim.calls))
	return mean(call -> call.responseTime, sim.calls) * (useMinutes ? 60*24 : 1) # if useMinutes, convert days to minutes
end

function getCallsReachedInTime(sim::Simulation;
	targetResponseTimes::Vector{Float} = sim.targetResponseTimes)
	@assert(sim.complete)
	@assert(all(call -> call.responseTime >= 0, sim.calls))
	return map(call -> call.responseTime <= targetResponseTimes[Int(call.priority)], sim.calls)
end

function countCallsReachedInTime(sim::Simulation;
	targetResponseTimes::Vector{Float} = sim.targetResponseTimes)
	@assert(sim.complete)
	@assert(all(call -> call.responseTime >= 0, sim.calls))
	return count(call -> call.responseTime <= targetResponseTimes[Int(call.priority)], sim.calls)
end

printDays(t::Float) = string(round(t, 2), " days")
printHours(t::Float) = string(round(t*24, 2), " hours") # convert days to hours, round
printMinutes(t::Float) = string(round(t*24*60, 2), " minutes") # convert days to minutes, round
printSeconds(t::Float) = string(round(t*24*60*60, 2), " seconds") # convert days to seconds, round
printPercent(p::Float) = string(round(p*100, 2), "%")

function printSimStats(sim::Simulation)
	println("=== Simulation statistics ===")
	println()
	println("Simulation duration = ", printDays(sim.endTime - sim.startTime))
	println()
	printAmbsStats(sim)
	println()
	printCallsStats(sim)
	println()
	printHospitalsStats(sim)
end

function printAmbsStats(sim::Simulation)
	@assert(sim.complete)
	
	ambsAvgDailyTravelTimes = [amb.totalTravelTime for amb in sim.ambulances] ./ (sim.time - sim.startTime)
	
	# Gadfly.plot(ecdf(ambsAvgDailyTravelTimes*24*60),
		# x="Average daily travel time (minutes)",
		# y="Cdf")
	
	println("Ambulance statistics:")
	
	println("Number of ambulances = ", sim.numAmbs)
	
	println("Average daily travel time: ")
	println(" mean = ", printHours(mean(ambsAvgDailyTravelTimes)))
	println(" std = ", printHours(std(ambsAvgDailyTravelTimes)))
	println(" min = ", printHours(minimum(ambsAvgDailyTravelTimes)))
	println(" max = ", printHours(maximum(ambsAvgDailyTravelTimes)))
end

function printCallsStats(sim::Simulation)
	@assert(sim.complete)
	
	responseTimes = getCallResponseTimes(sim)
	callsReachedInTime = getCallsReachedInTime(sim)
	
	println("Call statistics:")
	
	# # print, then remove, unanswered calls
	# callAnswered = [call.responseTime != nullTime for call in sim.calls] # may be false if sim not finished, or call cancelled
	# println("Unanswered calls: ", sim.numCalls - sum(callAnswered))
	# responseTimes = responseTimes[callAnswered] # no longer consider unanswered calls
	
	# Gadfly.plot(ecdf(responseTimes*24*60), x="Call response time (minutes)", y="Cdf")
	
	println("Number of calls = ", sim.numCalls)
	
	println("Response time: ")
	println(" mean = ", printMinutes(mean(responseTimes)))
	println(" std = ", printMinutes(std(responseTimes)))
	println(" min = ", printMinutes(minimum(responseTimes)))
	println(" max = ", printMinutes(maximum(responseTimes)))
	
	println("Reached in target response time: ")
	println(" mean = ", printPercent(mean(callsReachedInTime)))
	println(" sem = ", printPercent(sem(callsReachedInTime)))
end

function printHospitalsStats(sim::Simulation)
	@assert(sim.complete)
	
	println("Hospital statistics:")
	println("Number of hospitals = ", sim.numHospitals)
	println("Transfers:")
	hospitalsNumTransfers = map(h -> h.numTransfers, sim.hospitals)
	println(" mean = ", round(mean(hospitalsNumTransfers), 2))
	println(" std = ", round(std(hospitalsNumTransfers), 2))
end

"""
	function calcBatchMeans(values::Vector{Float}, batchSize::Int;
		batchGapSize::Int = 0, rmPartialBatch::Bool = false, returnBatchSizes::Bool = false)
Returns the batch means of `values` batched by size `batchSize`.

# Keyword arguments
- `batchGapSize` is the gap (number of values to ignore) between batches
- `rmPartialBatch` can be set to `true` to remove the last batch if it has size < `batchSize`
- `returnBatchSizes` can be set to `true` to also return the number of values in each batch
"""
function calcBatchMeans(values::Vector{Float}, batchSize::Int;
	batchGapSize::Int = 0, rmPartialBatch::Bool = false, returnBatchSizes::Bool = false)
	
	@assert(batchSize >= 1)
	@assert(batchGapSize >= 0)
	
	# calculate number of batches
	numValues = length(values)
	interBatchSize = batchSize + batchGapSize # number of values between start of subsequent batches
	numBatches = max(0, ceil(Int, numValues / interBatchSize))
	if rmPartialBatch
		numBatches = max(0, div(numValues + batchGapSize, interBatchSize))
	end
	
	# calculate batch means by time
	batchTotals = zeros(Float,numBatches) # total of values in each batch
	batchSizes = zeros(Int,numBatches) # total number of values in each batch
	for i = 1:numValues
		batchIndex = div(i-1, interBatchSize) + 1 # = ceil(Int, i / interBatchSize)
		isInGap = (i > batchIndex * interBatchSize - batchGapSize)
		if batchIndex >= 1 && batchIndex <= numBatches && !isInGap
			batchTotals[batchIndex] += values[i]
			batchSizes[batchIndex] += 1
		end
	end
	batchMeans = batchTotals ./ batchSizes
	
	# check batch sizes
	if rmPartialBatch
		@assert(all(n -> n == batchSize, batchSizes))
	else
		@assert(all(n -> n == batchSize, batchSizes[1:end-1]))
		@assert(numBatches == 0 || batchSizes[end] == rem(numValues, interBatchSize))
	end
	
	return returnBatchSizes ? (batchMeans, batchSizes) : batchMeans
end

"""
	function calcBatchMeans(values::Vector{Float}, times::Vector{Float}, batchTime::Float;
		startTime::Float = minimum(times), endTime::Float = maximum(times)*(1+eps(Float)),
		batchGapTime::Float = 0.0, rmPartialBatch::Bool = false, returnBatchSizes::Bool = false)
Returns the batch means of `values` batched by time, where `values[i]` corresponds with `times[i]`.
Each batch will batch values in a time interval `[t, t + batchTime)` for some `t`.
Empty batches have mean `NaN`.

# Keyword arguments
- `startTime` is the time at which batching starts; values with times before this are omitted
- `endTime` is the time at which batching ends; values with times at or after this are omitted
- `batchGapTime` is the time (duration) gap between batches
- `rmPartialBatch` can be set to `true` to remove the last batch if it does not end before `endTime - batchGapTime`.
- `returnBatchSizes` can be set to `true` to also return the number of values in each batch
"""
function calcBatchMeans(values::Vector{Float}, times::Vector{Float}, batchTime::Float;
	startTime::Float = minimum(times), endTime::Float = maximum(times)*(1+eps(Float)),
	batchGapTime::Float = 0.0, rmPartialBatch::Bool = false, returnBatchSizes::Bool = false)
	
	@assert(length(values) == length(times))
	@assert(batchTime > 0)
	@assert(startTime <= endTime)
	@assert(batchGapTime >= 0)
	
	# calculate number of batches
	interBatchTime = batchTime + batchGapTime # time difference between start of subsequent batches
	numBatches = max(0, ceil(Int, (endTime - startTime) / interBatchTime))
	if rmPartialBatch
		numBatches = max(0, floor(Int, (endTime - startTime + batchGapTime) / interBatchTime))
	end
	
	# calculate batch means by time
	batchTotals = zeros(Float,numBatches) # total of values in each batch
	batchSizes = zeros(Int,numBatches) # total number of values in each batch
	for (i,t) in enumerate(times)
		batchIndex = floor(Int, (t - startTime) / interBatchTime) + 1
		isInGap = (t - startTime >= batchIndex * interBatchTime - batchGapTime)
		if batchIndex >= 1 && batchIndex <= numBatches && !isInGap
			batchTotals[batchIndex] += values[i]
			batchSizes[batchIndex] += 1
		end
	end
	batchMeans = batchTotals ./ batchSizes
	
	return returnBatchSizes ? (batchMeans, batchSizes) : batchMeans
end

# calculate statistics on average response time based on batch means
function calcBatchMeanResponseTimes(sim::Simulation;
	batchTime::Float = nullTime, warmUpTime::Float = nullTime, coolDownTime::Float = nullTime,
	rmPartialBatch::Bool = true, returnBatchSizes::Bool = true)
	
	@assert(sim.complete == true)
	
	# collate data for use by calcBatchMeans function
	n = sim.numCalls
	times = [sim.calls[i].arrivalTime for i = 1:n]
	values = [sim.calls[i].responseTime for i = 1:n]
	
	return calcBatchMeans(values, times, batchTime;
		startTime = sim.startTime + warmUpTime, endTime = sim.endTime - coolDownTime,
		rmPartialBatch = rmPartialBatch, returnBatchSizes = returnBatchSizes)
end

# calculate Stats.sem for each row of x
function Stats.sem{T<:Real}(x::Array{T,2})
	return [Stats.sem(x[i,:]) for i = 1:size(x,1)]
end

# For each x[i], plot mean of y[i,:] with two sided confidence interval.
# Assumes that values in y[i,:] are from a normal distribution with unknown standard deviation.
function meanErrorPlot(x::Vector{Float}, y::Array{Float,2}, conf::Float=0.95)
	@assert(length(x) == size(y,1))
	@assert(0 < conf < 1)
	t = StatsFuns.tdistinvcdf(size(y,2)-1, 1-(1-conf)/2) # t-value, for two sided confidence interval
	# StatsFuns.tdistinvcdf(dof, p) # for n samples, dof = n - 1
	return Plots.scatter(x, mean(y,2), yerror=repmat(t*Stats.sem(y), 1, 2),
		xaxis=("x"), yaxis=("y"), m=(:hline), lab="")
end
function meanErrorPlot(x, y, conf::Float=0.95)
	meanErrorPlot(convert(Vector{Float},x), convert(Array{Float,2},y), conf)
end
# plot y values against indices (y[i,:] plotted at x = i)
function meanErrorPlot(y, conf::Float=0.95)
	meanErrorPlot(1:size(y,1), y, conf)
end

# For a vector of values, fit an AR(0) model, and return a p-value for the Durbin-Watson test.
# From wikipedia: "the Durbin–Watson statistic is a test statistic used to detect the presence
# of autocorrelation at lag 1 in the residuals (prediction errors) from a regression analysis".
function calcAR0DurbinWatsonTestPValue{T<:Real}(x::Vector{T})
	xFit = repmat([mean(x)], length(x), 1)
	residuals = x - xFit
	dwTest = HypothesisTests.DurbinWatsonTest(xFit, residuals)
	return HypothesisTests.pvalue(dwTest)
end
