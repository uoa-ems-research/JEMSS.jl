@testset "demand file i/o" begin
	demand = readDemandFile("data/demand/demand_1.csv")
	
	# check that demand can be written to file and read again
	tempDemandFilename = "temp/demand_1.csv"
	writeDemandFile(tempDemandFilename, demand)
	tempDemand = readDemandFile(tempDemandFilename)
	rm(tempDemandFilename)
	
	# check that input and output files have same raster filenames
	@test demand.rasterFilenames == tempDemand.rasterFilenames
	
	# check that demand files were read correctly
	for d in [demand, tempDemand]
		# check demand
		@test d.numRasters == 4
		@test d.numModes == 6
		@test d.numSets == 2
		
		# check demand modes
		@test [mode.index for mode in d.modes] == [1,2,3,4,5,6]
		@test [mode.rasterIndex for mode in d.modes] == [1,3,3,2,4,4]
		@test [Int(mode.priority) for mode in d.modes] == [1,2,3,1,2,3]
		@test [mode.arrivalRate for mode in d.modes] == [25,15,5,30,20,10]
		
		# check demand sets
		@test d.modeLookup == [1 2 3; 4 5 6]
		@test demand.setsStartTimes == sort(vcat(0,[0.7:1:10;],[0.3:1:10;]))
		@test demand.setsTimeOrder == [2-isodd(i) for i = 1:length(demand.setsTimeOrder)]
	end
end

@testset "demand mode lookup" begin
	demand = readDemandFile("data/demand/demand_1.csv")
	
	priorities = setdiff([instances(JEMSS.Priority)...], [JEMSS.nullPriority])
	for t = 0:0.1:1, priority in priorities
		i = JEMSS.getDemandMode!(demand, priority, t).index
		if 0.3 <= mod(t,1) < 0.7
			@test (priority == highPriority && i == 4) ||
				(priority == medPriority && i == 5) ||
				(priority == lowPriority && i == 6)
		else
			@test (priority == highPriority && i == 1) ||
				(priority == medPriority && i == 2) ||
				(priority == lowPriority && i == 3)
		end
	end
	
	# change in time beyond the last set start time should not cause demand set to change
	for t = demand.setsStartTimes[end] + (0:0.1:1)
		JEMSS.updateDemandToTime!(demand, t)
		@test demand.recentSetsStartTimesIndex == length(demand.setsStartTimes)
	end
end
