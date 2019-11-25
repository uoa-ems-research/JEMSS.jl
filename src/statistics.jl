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

function getCallResponseDurations(sim::Simulation)
	@assert(sim.complete)
	@assert(all(call -> call.responseDuration >= 0, sim.calls))
	return map(call -> call.responseDuration, sim.calls)
end

function getAvgCallResponseDuration(sim::Simulation; useMinutes::Bool = false)
	@assert(sim.complete)
	@assert(all(call -> call.responseDuration >= 0, sim.calls))
	return mean(call -> call.responseDuration, sim.calls) * (useMinutes ? 60*24 : 1) # if useMinutes, convert days to minutes
end

function getCallsReachedInTime(sim::Simulation;
	targetResponseDurations::Vector{Float} = sim.targetResponseDurations)
	@assert(sim.complete)
	@assert(all(call -> call.responseDuration >= 0, sim.calls))
	return map(call -> call.responseDuration <= targetResponseDurations[Int(call.priority)], sim.calls)
end

function countCallsReachedInTime(sim::Simulation;
	targetResponseDurations::Vector{Float} = sim.targetResponseDurations)
	@assert(sim.complete)
	@assert(all(call -> call.responseDuration >= 0, sim.calls))
	return count(call -> call.responseDuration <= targetResponseDurations[Int(call.priority)], sim.calls)
end

printDays(t::Float) = string(round(t, digits = 2), " days")
printHours(t::Float) = string(round(t*24, digits = 2), " hours") # convert days to hours, round
printMinutes(t::Float) = string(round(t*24*60, digits = 2), " minutes") # convert days to minutes, round
printSeconds(t::Float) = string(round(t*24*60*60, digits = 2), " seconds") # convert days to seconds, round
printPercent(p::Float) = string(round(p*100, digits = 2), "%")

printKilometres(d::Float) = string(round(d, digits = 2), " km")

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
	
	ambsAvgDailyTravelTimes = [amb.statusDurations[ambTravelling] for amb in sim.ambulances] ./ (sim.time - sim.startTime)
	ambsAvgDailyTravelDists = [amb.statusDistances[ambTravelling] for amb in sim.ambulances] ./ (sim.time - sim.startTime)
	
	# Gadfly.plot(ecdf(ambsAvgDailyTravelTimes*24*60),
		# x="Average daily travel time (minutes)",
		# y="Cdf")
	
	println("Ambulance statistics:")
	
	println("Number of ambulances = ", sim.numAmbs)
	
	println("Average daily travel time: ")
	println(" mean = ", printHours(mean(ambsAvgDailyTravelTimes)))
	println(" std = ", printHours(std(ambsAvgDailyTravelTimes, corrected = false)))
	println(" min = ", printHours(minimum(ambsAvgDailyTravelTimes)))
	println(" max = ", printHours(maximum(ambsAvgDailyTravelTimes)))
	
	println("Average daily travel distance: ")
	println(" mean = ", printKilometres(mean(ambsAvgDailyTravelDists)))
	println(" std = ", printKilometres(std(ambsAvgDailyTravelDists, corrected = false)))
	println(" min = ", printKilometres(minimum(ambsAvgDailyTravelDists)))
	println(" max = ", printKilometres(maximum(ambsAvgDailyTravelDists)))
end

function printCallsStats(sim::Simulation)
	@assert(sim.complete)
	
	responseDurations = getCallResponseDurations(sim)
	callsReachedInTime = getCallsReachedInTime(sim)
	
	println("Call statistics:")
	
	# # print, then remove, unanswered calls
	# callAnswered = [call.responseDuration != nullTime for call in sim.calls] # may be false if sim not finished, or call cancelled
	# println("Unanswered calls: ", sim.numCalls - sum(callAnswered))
	# responseDurations = responseDurations[callAnswered] # no longer consider unanswered calls
	
	# Gadfly.plot(ecdf(responseDurations*24*60), x="Call response duration (minutes)", y="Cdf")
	
	println("Number of calls = ", sim.numCalls)
	
	println("Response duration: ")
	println(" mean = ", printMinutes(mean(responseDurations)))
	println(" std = ", printMinutes(std(responseDurations, corrected = false)))
	println(" min = ", printMinutes(minimum(responseDurations)))
	println(" max = ", printMinutes(maximum(responseDurations)))
	
	println("Reached in target response duration: ")
	println(" mean = ", printPercent(mean(callsReachedInTime)))
	println(" sem = ", printPercent(sem(callsReachedInTime)))
end

function printHospitalsStats(sim::Simulation)
	@assert(sim.complete)
	
	println("Hospital statistics:")
	println("Number of hospitals = ", sim.numHospitals)
	println("Calls:")
	hospitalsNumCalls = map(h -> h.numCalls, sim.hospitals)
	println(" mean = ", round(mean(hospitalsNumCalls), digits = 2))
	println(" std = ", round(std(hospitalsNumCalls, corrected = false), digits = 2))
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

# calculate statistics on average response duration based on batch means
function calcBatchMeanResponseDurations(sim::Simulation;
	batchTime::Float = nullTime, warmUpTime::Float = nullTime, coolDownTime::Float = nullTime,
	rmPartialBatch::Bool = true, returnBatchSizes::Bool = true)
	
	@assert(sim.complete == true)
	
	# collate data for use by calcBatchMeans function
	n = sim.numCalls
	times = [sim.calls[i].arrivalTime for i = 1:n]
	values = [sim.calls[i].responseDuration for i = 1:n]
	
	return calcBatchMeans(values, times, batchTime;
		startTime = sim.startTime + warmUpTime, endTime = sim.endTime - coolDownTime,
		rmPartialBatch = rmPartialBatch, returnBatchSizes = returnBatchSizes)
end

# calculate sem for 2d array, calculating along given dim
function semDim(x::Array{T,2}, dim::Int = 1) where T <: Real
	@assert(dim == 1 || dim == 2)
	dim == 1 && return collect([sem(x[:,j]) for j = 1:size(x,2)]') # sem of each col
	dim == 2 && return [sem(x[i,:]) for i = 1:size(x,1)] # sem of each row
end

# For each x[i], plot mean of y[i,:] with two sided confidence interval.
# Assumes that values in y[i,:] are from a normal distribution with unknown standard deviation.
function meanErrorPlot(x::Vector{Float}, y::Array{Float,2}, conf::Float=0.95)
	@assert(length(x) == size(y,1))
	@assert(0 < conf < 1)
	t = StatsFuns.tdistinvcdf(size(y,2)-1, 1-(1-conf)/2) # t-value, for two sided confidence interval
	# StatsFuns.tdistinvcdf(dof, p) # for n samples, dof = n - 1
	@warn("meanErrorPlot still requires use of Plots package; need to change this to RecipesBase")
	return Plots.scatter(x, mean(y, dims = 2), yerror=repeat(t*semDim(y,2), 1, 2),
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
function calcAR0DurbinWatsonTestPValue(x::Vector{T}) where T <: Real
	xFit = repeat([mean(x)], length(x))
	residuals = x - xFit
	dwTest = HypothesisTests.DurbinWatsonTest(xFit, residuals)
	return HypothesisTests.pvalue(dwTest)
end

# Return half-width of two-sided confidence interval of estimate of mean.
# Samples should come from a normal distribution, standard deviation is assumed to be unknown.
function tDistrHalfWidth(x::Vector{T}; conf = 0.95) where T <: Real
	t = StatsFuns.tdistinvcdf(length(x)-1, 1-(1-conf)/2) # t-value, for two sided confidence interval
	return t * sem(x)
end

# confidence interval
function confInterval(mhw::MeanAndHalfWidth)
	return (mhw.mean - mhw.halfWidth, mhw.mean + mhw.halfWidth)
end

# Return a list of the main stats periods.
# Assumes that the first period (after any warm-up) is the main period.
# Can only return multiple periods if `sim.stats.periodDurationsIter` is cyclical.
# Useful for batches in a single simulation replication.
function getPeriodStatsList(sim::Simulation)::Vector{SimPeriodStats}
	@assert(sim.complete)
	stats = sim.stats # shorthand
	i = stats.warmUpDuration == 0 ? 1 : 2 # ignore first period if it is for warm-up
	itr = stats.periodDurationsIter.itr # can be empty if stats control not set
	k = isa(itr, Iterators.Cycle) ? length(itr.xs) : length(itr) # gap between period indices
	periods = stats.periods[i:k:end]
	if isempty(periods) return periods end
	periodDuration = periods[1].duration
	if !isapprox(periods[end].duration, periodDuration) pop!(periods) end
	@assert(all(p -> isapprox(p.duration, periodDuration), periods))
	return periods
end

# given the sim replications, return a list with the main stats period for each rep
function getRepsPeriodStatsList(reps::Vector{Simulation}; periodIndex::Int = nullIndex)::Vector{SimPeriodStats}
	@assert(all(rep -> rep.complete, reps))
	if periodIndex == nullIndex
		warmUpDuration = reps[1].stats.warmUpDuration
		@assert(all(rep -> rep.stats.warmUpDuration == warmUpDuration, reps))
		periodIndex = warmUpDuration == 0 ? 1 : 2
	end
	periods = [rep.stats.periods[periodIndex] for rep in reps]
	@assert(all(p -> isapprox(p.duration, periods[1].duration), periods))
	return periods
end

# given a sim replication, return the main stats period
function getRepPeriodStats(rep::Simulation; periodIndex::Int = nullIndex)::SimPeriodStats
	return getRepsPeriodStatsList([rep], periodIndex = periodIndex)[1]
end

# Return a dictionary of statistics from a list of period statistics.
# Periods should be the same duration.
# Confidence intervals assume that samples obtained from periods are IID, and are from a population with a normal distribution and unknown standard deviation.
# The function flatten(dict) may be useful if writing this dict to file.
function statsDictFromPeriodStatsList(periods::Vector{SimPeriodStats}; conf = 0.95)
	
	if isempty(periods)
		period = SimPeriodStats()
		period.duration = 0.0
		ambStatuses = (instances(AmbStatus)..., instances(AmbStatusSet)...)
		period.ambulance.statusDurations = Dict([s => 0.0 for s in ambStatuses])
		period.ambulance.statusDistances = Dict([s => 0.0 for s in ambStatuses])
		period.callPriorities = Dict([p => CallStats() for p in priorities])
		periods = [period]
	end
	
	duration = periods[1].duration
	ambDays = length(periods[1].ambulances) * duration
	@assert(all(p -> isapprox(p.duration, duration), periods))
	
	d = statsDict = Dict{String,Any}()
	
	meanAndHalfWidth(x; conf = conf) = MeanAndHalfWidth(mean(x), tDistrHalfWidth(x; conf = conf))
	
	###########
	# ambulance
	
	d["ambs"] = Dict{String,Any}()
	
	function getAmbsStat(statName::Symbol)
		x = [getfield(p.ambulance, statName) for p in periods] / ambDays
		return meanAndHalfWidth(x)
	end
	function getAmbsDurationStat(s::Union{AmbStatus,AmbStatusSet})
		x = [p.ambulance.statusDurations[s] for p in periods] / ambDays
		return meanAndHalfWidth(x)
	end
	function getAmbsDistanceStat(s::Union{AmbStatus,AmbStatusSet})
		x = [p.ambulance.statusDistances[s] for p in periods] / ambDays
		return meanAndHalfWidth(x)
	end
	
	# counts
	d["ambs"]["avgDailyNumCallsTreated"] = getAmbsStat(:numCallsTreated)
	d["ambs"]["avgDailyNumCallsTransported"] = getAmbsStat(:numCallsTransported)
	d["ambs"]["avgDailyNumDispatches"] = getAmbsStat(:numDispatches)
	d["ambs"]["avgDailyNumDispatchesFromStation"] = getAmbsStat(:numDispatchesFromStation)
	d["ambs"]["avgDailyNumDispatchesOnRoad"] = getAmbsStat(:numDispatchesOnRoad)
	d["ambs"]["avgDailyNumDispatchesOnFree"] = getAmbsStat(:numDispatchesOnFree)
	d["ambs"]["avgDailyNumRedispatches"] = getAmbsStat(:numRedispatches)
	d["ambs"]["avgDailyNumMoveUps"] = getAmbsStat(:numMoveUps)
	d["ambs"]["avgDailyNumMoveUpsFromStation"] = getAmbsStat(:numMoveUpsFromStation)
	d["ambs"]["avgDailyNumMoveUpsOnRoad"] = getAmbsStat(:numMoveUpsOnRoad)
	d["ambs"]["avgDailyNumMoveUpsOnFree"] = getAmbsStat(:numMoveUpsOnFree)
	d["ambs"]["avgDailyNumMoveUpsReturnToPrevStation"] = getAmbsStat(:numMoveUpsReturnToPrevStation)
	
	# durations - status
	d["ambs"]["avgDailySleepingDurationHours"] = getAmbsDurationStat(ambSleeping) * 24
	d["ambs"]["avgDailyIdleAtStationDurationHours"] = getAmbsDurationStat(ambIdleAtStation) * 24
	d["ambs"]["avgDailyGoingToCallDurationHours"] = getAmbsDurationStat(ambGoingToCall) * 24
	d["ambs"]["avgDailyAtCallDurationHours"] = getAmbsDurationStat(ambAtCall) * 24
	d["ambs"]["avgDailyGoingToHospitalDurationHours"] = getAmbsDurationStat(ambGoingToHospital) * 24
	d["ambs"]["avgDailyAtHospitalDurationHours"] = getAmbsDurationStat(ambAtHospital) * 24
	d["ambs"]["avgDailyReturningToStationDurationHours"] = getAmbsDurationStat(ambReturningToStation) * 24
	d["ambs"]["avgDailyMovingUpToStationDurationHours"] = getAmbsDurationStat(ambMovingUpToStation) * 24
	
	# durations - status sets
	d["ambs"]["avgDailyWorkingDurationHours"] = getAmbsDurationStat(ambWorking) * 24
	d["ambs"]["avgDailyBusyDurationHours"] = getAmbsDurationStat(ambBusy) * 24
	d["ambs"]["avgDailyFreeDurationHours"] = getAmbsDurationStat(ambFree) * 24
	d["ambs"]["avgDailyTravelDurationHours"] = getAmbsDurationStat(ambTravelling) * 24
	d["ambs"]["avgDailyGoingToStationDurationHours"] = getAmbsDurationStat(ambGoingToStation) * 24
	
	# distances
	d["ambs"]["avgDailyTravelDistanceKms"] = getAmbsDistanceStat(ambTravelling)
	d["ambs"]["avgDailyGoingToCallDistanceKms"] = getAmbsDistanceStat(ambGoingToCall)
	d["ambs"]["avgDailyGoingToHospitalDistanceKms"] = getAmbsDistanceStat(ambGoingToHospital)
	d["ambs"]["avgDailyReturningToStationDistanceKms"] = getAmbsDistanceStat(ambReturningToStation)
	d["ambs"]["avgDailyMovingUpToStationDistanceKms"] = getAmbsDistanceStat(ambMovingUpToStation)
	
	###########
	# call
	
	d["calls"] = Dict{String,Any}()
	
	for priority in instances(Priority) # will use nullPriority to indicate all priorities
		name = priority == nullPriority ? "all" : string(priority)
		d1 = d["calls"][name] = Dict{String,Any}()
		
		callStatsList = priority == nullPriority ? [p.call for p in periods] : [p.callPriorities[priority] for p in periods]
		function getCallsStat(dividendName::Symbol, divisorName::Symbol = :numCalls)
			x = [getfield(callStats, dividendName) for callStats in callStatsList]
			y = [getfield(callStats, divisorName) for callStats in callStatsList]
			i = findall(x -> x != 0, y)
			# if isempty(i) return MeanAndHalfWidth(0,0) end # do not do this; lack of samples doesn't imply mean of zero
			return meanAndHalfWidth(x[i]./y[i])
		end
		
		# fractions (between 0 and 1); don't need to write "avgFrac", as fraction is already an average
		d1["fracQueued"] = getCallsStat(:numQueued)
		d1["fracBumped"] = getCallsStat(:numBumped)
		d1["fracTransported"] = getCallsStat(:numTransports)
		d1["fracResponsesInTime"] = getCallsStat(:numResponsesInTime)
		
		# durations
		d1["avgDispatchDelayMinutes"] = getCallsStat(:totalDispatchDelay) * (24*60)
		d1["avgOnSceneDurationMinutes"] = getCallsStat(:totalOnSceneDuration) * (24*60)
		d1["avgHandoverDurationMinutes"] = getCallsStat(:totalHandoverDuration, :numTransports) * (24*60)
		d1["avgQueuedDurationMinutes"] = getCallsStat(:totalQueuedDuration, :numQueued) * (24*60)
		d1["avgBumpedDurationMinutes"] = getCallsStat(:totalBumpedDuration, :numBumped) * (24*60)
		d1["avgWaitingForAmbDurationMinutes"] = getCallsStat(:totalWaitingForAmbDuration) * (24*60)
		d1["avgResponseDurationMinutes"] = getCallsStat(:totalResponseDuration) * (24*60)
		d1["avgAmbGoingToCallDurationMinutes"] = getCallsStat(:totalAmbGoingToCallDuration) * (24*60)
		d1["avgTransportDurationMinutes"] = getCallsStat(:totalTransportDuration, :numTransports) * (24*60)
		d1["avgServiceDurationMinutes"] = getCallsStat(:totalServiceDuration) * (24*60)
		
		# misc
		d1["avgNumBumps"] = getCallsStat(:numBumps)
		d1["avgDailyNumCalls"] = meanAndHalfWidth([callStats.numCalls / duration for callStats in callStatsList])
	end
	
	###########
	
	return statsDict
end
