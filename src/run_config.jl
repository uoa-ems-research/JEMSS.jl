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

# run configuration xml file
function runConfig(configFilename::String; doPrint::Bool=false)
    sim = initSim(configFilename; allowResim=true, allowWriteOutput=true, doPrint=doPrint)
    if !isempty(sim.reps)
        simulateReps!(sim)
    else
        openOutputFiles!(sim)
        simulate!(sim)
        closeOutputFiles!(sim)
        doPrint && printSimStats(sim)
    end
    if sim.writeOutput
        writeOutputFiles(sim)
        sim.writeOutput = sim.backup.writeOutput = false # need to switch off in order to re-run the simulation
    end
    return sim
end

# initialise simulation from input files
function initSim(configFilename::String;
    allowResim::Bool=false, createBackup::Bool=true, allowWriteOutput::Bool=false, doPrint::Bool=false)

    # read sim config xml file
    configFilename = configFilename |> interpolateString |> abspath
    global configFileDir = dirname(configFilename)
    rootElt = xmlFileRoot(configFilename)
    @assert(xName(rootElt) == "simConfig", string("xml root has incorrect name: ", xName(rootElt)))

    # for progress messages:
    t = 0.0
    printInitMsg(msg) = doPrint && (t = time(); print(msg))
    printInitTime() = doPrint && println(": ", round(time() - t, digits=2), " seconds")

    ##################
    # sim config

    printInitMsg("reading config file data")

    sim = Simulation()
    sim.configRootElt = rootElt

    # input
    sim.inputPath = joinPathIfNotAbs(configFileDir, eltContentInterpVal(rootElt, "inputPath")) # input path can be absolute, or relative to configFilename
    simFilesElt = findElt(rootElt, "simFiles")
    inputFiles = childrenNodeNames(simFilesElt)
    sim.inputFiles = Dict{String,File}()
    for inputFile in inputFiles
        file = File()
        file.path = joinPathIfNotAbs(sim.inputPath, eltContentInterpVal(simFilesElt, inputFile))
        file.name = basename(file.path)
        if inputFile != "rNetTravels" # do not need checksum of rNetTravels file
            file.checksum = fileChecksum(file.path)
        end
        sim.inputFiles[inputFile] = file
    end

    # output
    sim.writeOutput = allowWriteOutput && eltContentVal(rootElt, "writeOutput")
    sim.outputPath = joinPathIfNotAbs(configFileDir, eltContentInterpVal(rootElt, "outputPath")) # output path can be absolute, or relative to configFilename
    isdir(sim.outputPath) || mkdir(sim.outputPath)
    outputFilesElt = findElt(rootElt, "outputFiles")
    outputFiles = childrenNodeNames(outputFilesElt)
    sim.outputFiles = Dict{String,File}()
    for outputFile in outputFiles
        file = File()
        file.path = joinPathIfNotAbs(sim.outputPath, eltContentInterpVal(outputFilesElt, outputFile))
        file.name = basename(file.path)
        sim.outputFiles[outputFile] = file
    end

    # events file filter
    if containsElt(rootElt, "eventsFileFilter")
        eventsFiltered = eltContentVal(rootElt, "eventsFileFilter")
        @assert(isa(eventsFiltered, Vector{EventForm}))
        for eventForm in setdiff(instances(EventForm), eventsFiltered)
            sim.eventsFile.eventFilter[eventForm] = false
        end
    end

    printInitTime()

    ##################
    # read simulation input files

    printInitMsg("reading input data")

    simFilePath(name::String) = sim.inputFiles[name].path
    hasInputFile(name::String) = haskey(sim.inputFiles, name)

    # read sim data
    sim.ambulances = readAmbsFile(simFilePath("ambulances"))
    sim.hospitals = readHospitalsFile(simFilePath("hospitals"))
    sim.stations = readStationsFile(simFilePath("stations"))

    # calls
    callGenConfig = nothing
    @assert(hasInputFile("calls") + hasInputFile("callGenConfig") == 1, "Need exactly one of these input files: calls, callGenConfig.")
    if hasInputFile("calls")
        (sim.calls, sim.startTime) = readCallsFile(simFilePath("calls"))
    elseif hasInputFile("callGenConfig")
        callGenConfig = readGenConfig(simFilePath("callGenConfig"))
        sim.startTime = callGenConfig.startTime
    else
        error()
    end
    sim.time = sim.startTime

    sim.numAmbs = length(sim.ambulances)
    sim.numCalls = length(sim.calls)
    sim.numHospitals = length(sim.hospitals)
    sim.numStations = length(sim.stations)

    # read network data
    sim.net = Network()
    net = sim.net # shorthand
    fGraph = net.fGraph # shorthand
    fGraph.nodes = readNodesFile(simFilePath("nodes"))
    (fGraph.arcs, arcTravelTimes) = readArcsFile(simFilePath("arcs"))

    # read rNetTravels from file, if saved
    rNetTravelsLoaded = NetTravel[]
    rNetTravelsFilename = ""
    if hasInputFile("rNetTravels")
        rNetTravelsFilename = simFilePath("rNetTravels")
        if isfile(rNetTravelsFilename)
            try
                rNetTravelsLoaded = readRNetTravelsFile(rNetTravelsFilename)
            catch
                @warn("failed to read rNetTravels file")
            end
        elseif !isdir(dirname(rNetTravelsFilename)) || basename(rNetTravelsFilename) == ""
            # rNetTravelsFilename is invalid
            rNetTravelsFilename = ""
        else
            # save net.rNetTravels to file once calculated
        end
    end

    # read misc
    sim.map = readMapFile(simFilePath("map"))
    map = sim.map # shorthand
    (sim.targetResponseDurations, sim.responseTravelPriorities) = readPrioritiesFile(simFilePath("priorities"))
    sim.travel = readTravelFile(simFilePath("travel"))
    # "demand" file can be slow to read, will not read here but read elsewhere when needed
    if hasInputFile("demandCoverage")
        sim.demandCoverage = readDemandCoverageFile(simFilePath("demandCoverage"))
    end
    if hasInputFile("mobilisationDelay")
        sim.mobilisationDelay = readMobilisationDelayFile(simFilePath("mobilisationDelay"))
    end

    printInitTime()

    ##################
    # grid

    printInitMsg("placing nodes in grid")

    # hard-coded grid size
    # grid rects will be roughly square, with one node per square on average
    n = length(fGraph.nodes)
    xDist = map.xRange * map.xScale
    yDist = map.yRange * map.yScale
    nx = Int(ceil(sqrt(n * xDist / yDist)))
    ny = Int(ceil(sqrt(n * yDist / xDist)))

    sim.grid = Grid(map, nx, ny)
    grid = sim.grid # shorthand
    gridPlaceNodes!(map, grid, fGraph.nodes)
    printInitTime()

    doPrint && println("nodes: ", length(fGraph.nodes), ", grid size: ", nx, " x ", ny)

    ##################
    # generate calls

    callSets = [Call[]] # init
    if hasInputFile("callGenConfig")
        printInitMsg("generating calls")

        numReps = sim.numReps = callGenConfig.numCallsFiles # shorthand
        @assert(numReps >= 1)
        callSets = [makeCalls(callGenConfig) for i = 1:numReps]

        if numReps > 1
            setSimReps!(sim, callSets)
        end

        sim.calls = callSets[1] # set sim.calls for first replication
        sim.numCalls = length(sim.calls)

        printInitTime()

        (numReps > 1) && doPrint && println("number of call sets: $numReps")
    end

    ##################
    # network

    printInitMsg("initialising fGraph")
    initGraph!(fGraph)
    printInitTime()

    if any(arc -> isnan(arc.distance), fGraph.arcs)
        printInitMsg("calculating arc distances for arcs with NaN distance")
        setArcDistances!(fGraph, map)
        printInitTime()
    end

    printInitMsg("checking fGraph")
    checkGraph(fGraph, map)
    printInitTime()

    printInitMsg("initialising fNetTravels")
    initFNetTravels!(net, arcTravelTimes)
    printInitTime()

    printInitMsg("creating rGraph from fGraph")
    createRGraphFromFGraph!(net)
    printInitTime()
    doPrint && println("fNodes: ", length(net.fGraph.nodes), ", rNodes: ", length(net.rGraph.nodes))

    printInitMsg("checking rGraph")
    checkGraph(net.rGraph, map)
    printInitTime()

    if rNetTravelsLoaded != []
        doPrint && println("using data from rNetTravels file")
        try
            printInitMsg("creating rNetTravels from fNetTravels")
            createRNetTravelsFromFNetTravels!(net; rNetTravelsLoaded=rNetTravelsLoaded)
            printInitTime()
        catch
            doPrint && println()
            @warn("failed to use data from rNetTravels file")
            rNetTravelsLoaded = []
        end
    end
    if rNetTravelsLoaded == []
        printInitMsg("creating rNetTravels from fNetTravels, and shortest paths")
        createRNetTravelsFromFNetTravels!(net)
        printInitTime()
        if rNetTravelsFilename != ""
            printInitMsg("saving rNetTravels to file")
            writeRNetTravelsFile(rNetTravelsFilename, net.rNetTravels)
            printInitTime()
        end
    end

    ##################
    # travel

    printInitMsg("initialising travel")

    travel = sim.travel # shorthand
    @assert(travel.setsStartTimes[1] <= sim.startTime)
    @assert(length(net.fNetTravels) == travel.numModes)
    for travelMode in travel.modes
        travelMode.fNetTravel = net.fNetTravels[travelMode.index]
        travelMode.rNetTravel = net.rNetTravels[travelMode.index]
    end

    printInitTime()

    ##################
    # sim - ambulances, calls, hospitals, stations...

    printInitMsg("adding ambulances, calls, etc")

    # for each call, hospital, and station, find nearest node
    for calls in (sim.calls, callSets[2:end]...) # note that sim.calls should be callSets[1] (if callSets is not empty)
        for c in calls
            if c.nearestNodeIndex == nullIndex
                (c.nearestNodeIndex, c.nearestNodeDist) = findNearestNode(map, grid, fGraph.nodes, c.location)
            end
        end
    end
    for h in sim.hospitals
        (h.nearestNodeIndex, h.nearestNodeDist) = findNearestNode(map, grid, fGraph.nodes, h.location)
    end
    for s in sim.stations
        (s.nearestNodeIndex, s.nearestNodeDist) = findNearestNode(map, grid, fGraph.nodes, s.location)
    end

    # create event list
    # try to add events to eventList in reverse time order, to reduce sorting required
    sim.eventList = Vector{Event}()

    # add first call to event list
    addEvent!(sim.eventList, sim.calls[1])

    # create ambulance wake up events
    for a in sim.ambulances
        initAmbulance!(sim, a)
        # currently, this sets ambulances to wake up at start of sim, since wake up and sleep events are not in ambulances file yet
    end

    # init station stats
    for station in sim.stations
        station.numIdleAmbsTotalDuration = OffsetVector(zeros(Float, sim.numAmbs + 1), 0:sim.numAmbs)
        station.currentNumIdleAmbsSetTime = sim.startTime
    end

    printInitTime()

    printInitMsg("storing times between fNodes and common locations")

    # for each station, find time to each node in fGraph for each travel mode, and vice versa (node to station)
    # for each node in fGraph and each travel mode, find nearest hospital
    # requires deterministic and static travel times

    commonFNodes = sort(unique(vcat([h.nearestNodeIndex for h in sim.hospitals], [s.nearestNodeIndex for s in sim.stations])))
    setCommonFNodes!(net, commonFNodes)

    # find the nearest hospital to travel to from each node in fGraph
    numFNodes = length(fGraph.nodes) # shorthand
    for fNetTravel in net.fNetTravels
        fNetTravel.fNodeNearestHospitalIndex = Vector{Int}(undef, numFNodes)
        travelModeIndex = fNetTravel.modeIndex # shorthand
        travelMode = travel.modes[travelModeIndex] # shorthand
        for node in fGraph.nodes
            # find nearest hospital to node
            minTime = Inf
            nearestHospitalIndex = nullIndex
            for hospital in sim.hospitals
                travelTime = shortestPathTravelTime(net, travelModeIndex, node.index, hospital.nearestNodeIndex)
                travelTime += offRoadTravelTime(travelMode, hospital.nearestNodeDist)
                if travelTime < minTime
                    minTime = travelTime
                    nearestHospitalIndex = hospital.index
                end
            end
            fNetTravel.fNodeNearestHospitalIndex[node.index] = nearestHospitalIndex
        end
    end

    printInitTime()

    ##################
    # decision logic

    decisionElt = findElt(rootElt, "decision")
    sim.addCallToQueue! = eltContentVal(decisionElt, "callQueueing")
    sim.findAmbToDispatch! = eltContentVal(decisionElt, "dispatch")
    hasInputFile("redispatch") && (sim.redispatch = readRedispatchFile(simFilePath("redispatch")))

    # move up
    mud = sim.moveUpData # shorthand
    moveUpElt = findElt(decisionElt, "moveUp")
    moveUpModuleName = eltContent(moveUpElt, "module")

    if moveUpModuleName == "none"
        mud.useMoveUp = false
        mud.moveUpModule = nullMoveUpModule
        doPrint && println("not using move up")
    else
        printInitMsg("initialising move up")

        mud.useMoveUp = true

        if moveUpModuleName == "comp_table"
            mud.moveUpModule = compTableModule
            compTableElt = findElt(moveUpElt, "compTable")
            compTableFilename = joinPathIfNotAbs(sim.inputPath, eltContentInterpVal(compTableElt, "filename"))
            initCompTable!(sim, compTableFilename)

        elseif moveUpModuleName == "ddsm"
            mud.moveUpModule = ddsmModule
            ddsmElt = findElt(moveUpElt, "ddsm")
            initDdsm!(sim;
                coverFractionTargetT1=eltContentVal(ddsmElt, "coverFractionTargetT1"),
                travelTimeCost=eltContentVal(ddsmElt, "travelTimeCost"),
                slackWeight=eltContentVal(ddsmElt, "slackWeight"),
                coverTimeDemandPriorities=eltContentVal(ddsmElt, "coverTimeDemandPriorities"),
                options=convert(Dict{Symbol,Any}, Dict(eltContentVal(ddsmElt, "options"))))

        elseif moveUpModuleName == "dmexclp"
            mud.moveUpModule = dmexclpModule
            dmexclpElt = findElt(moveUpElt, "dmexclp")
            initDmexclp!(sim; busyFraction=eltContentVal(dmexclpElt, "busyFraction"))

        elseif moveUpModuleName == "multi_comp_table"
            mud.moveUpModule = multiCompTableModule
            multiCompTableElt = findElt(moveUpElt, "multiCompTable")
            multiCompTableFilename = joinPathIfNotAbs(sim.inputPath, eltContentInterpVal(multiCompTableElt, "filename"))
            initMultiCompTable!(sim, multiCompTableFilename)

        elseif moveUpModuleName == "priority_list"
            mud.moveUpModule = priorityListModule
            priorityListElt = findElt(moveUpElt, "priorityList")
            priorityListFilename = joinPathIfNotAbs(sim.inputPath, eltContentInterpVal(priorityListElt, "filename"))
            initPriorityList!(sim, priorityListFilename)

        elseif moveUpModuleName == "zhang_ip"
            mud.moveUpModule = zhangIpModule
            zhangIpElt = findElt(moveUpElt, "zhangIp")
            zhangIpParamsFilename = joinPathIfNotAbs(sim.inputPath, eltContentInterpVal(zhangIpElt, "paramsFilename"))
            initZhangIp!(sim;
                paramsFilename=zhangIpParamsFilename)

        elseif moveUpModuleName == "temp0"
            mud.moveUpModule = temp0Module
            temp0Elt = findElt(moveUpElt, "temp0")
            initTemp0!(sim;
                busyFraction=eltContentVal(temp0Elt, "busyFraction"),
                travelTimeCost=eltContentVal(temp0Elt, "travelTimeCost"),
                maxIdleAmbTravelTime=eltContentVal(temp0Elt, "maxIdleAmbTravelTime"),
                maxNumNearestStations=eltContentVal(temp0Elt, "maxNumNearestStations"))

        elseif moveUpModuleName == "temp1"
            mud.moveUpModule = temp1Module
            temp1Elt = findElt(moveUpElt, "temp1")
            initTemp1!(sim;
                busyFraction=eltContentVal(temp1Elt, "busyFraction"),
                travelTimeCost=eltContentVal(temp1Elt, "travelTimeCost"),
                maxIdleAmbTravelTime=eltContentVal(temp1Elt, "maxIdleAmbTravelTime"),
                maxNumNearestStations=eltContentVal(temp1Elt, "maxNumNearestStations"))

        elseif moveUpModuleName == "temp2"
            mud.moveUpModule = temp2Module
            temp2Elt = findElt(moveUpElt, "temp2")
            initTemp2!(sim;
                busyFraction=eltContentVal(temp2Elt, "busyFraction"),
                travelTimeCost=eltContentVal(temp2Elt, "travelTimeCost"),
                maxIdleAmbTravelTime=eltContentVal(temp2Elt, "maxIdleAmbTravelTime"),
                maxNumNearestStations=eltContentVal(temp2Elt, "maxNumNearestStations"))
        else
            error("invalid move up module name given: ", moveUpModuleName)
        end

        printInitTime()

        doPrint && println("using move up module: ", moveUpModuleName)
    end

    ##################
    # resimulation

    if allowResim
        sim.resim.use = eltContentVal(rootElt, "resim")
        if sim.resim.use
            printInitMsg("")
            initResim!(sim, doPrint=doPrint)
            doPrint && print("initialised resimulation")
            printInitTime()
        end
    end

    ##################
    # statistics

    if hasInputFile("statsControl")
        printInitMsg("initialising statistics")
        stats = sim.stats # shorthand
        stats.doCapture = true
        dict = readStatsControlFile(simFilePath("statsControl"))
        for fname in (:periodDurationsIter, :warmUpDuration, :recordResponseDurationHist, :responseDurationHistBinWidth)
            setfield!(stats, fname, dict[string(fname)])
        end
        stats.nextCaptureTime = sim.startTime + (stats.warmUpDuration > 0 ? stats.warmUpDuration : first(stats.periodDurationsIter))
        printInitTime()
    end

    ##################
    # misc

    sim.initialised = true # at this point, the simulation could be run

    # set sim replications
    if sim.numReps >= 1
        makeRepsRunnable!(sim) # need to do this last in sim init (but can do before backing up), as it copies fields from sim
    end

    ##################
    # sim backup

    if createBackup
        printInitMsg("creating sim backup")
        backup!(sim) # for restarting sim
        printInitTime()
    end

    return sim
end
