__precompile__()
module JEMSS

# animation
using HttpServer
using WebSockets
using JSON

# files
using LightXML

# optimisation (move-up)
using JuMP
using GLPKMathProgInterface # does not use precompile

# statistics
using Distributions
using HypothesisTests
using Stats
using StatsFuns

# misc
using LightGraphs
using ArchGDAL # does not use precompile
import Plots

# simulation functions
export
	initSimulation, runConfig, # run_config
	simulate!, simulateToTime!, simulateToEnd!, backupSim!, resetSim!, simulateNextEvent!, # simulation
	animate # animation

# file functions
export
	readDlmFileNextLine!, readDlmFile, openNewFile, writeDlmLine!, arrayDict, writeTablesToFile!, writeTablesToFile, readTablesFromFile, readTablesFromData, tableRowsFieldDicts, fileChecksum, serializeToFile, deserializeFile, interpolateString, xmlFileRoot, findElt, eltContent, eltContentVal, eltContentInterpVal, childrenNodeNames, selectXmlFile, # file_io
	runGenConfig, # gen_sim_files
	readAmbsFile, readArcsFile, readCallsFile, readCompTableFile, readDemandFile, readEventsFile, readGeoFile, readHospitalsFile, readMapFile, readNodesFile, readPrioritiesFile, readPriorityListFile, readRasterFile, readRNetTravelsFile, readStationsFile, readTravelFile, readDeploymentPoliciesFile, readZhangIpParamsFile, # read_sim_files
	writeAmbsFile, writeArcsFile, writeCallsFile, writeDemandFile, writeHospitalsFile, writeMapFile, writeNodesFile, writePrioritiesFile, writeRNetTravelsFile, writeStationsFile, writeTravelFile, openOutputFiles!, closeOutputFiles!, writeEventToFile!, writeStatsFiles!, writeDeploymentPoliciesFile, writeBatchMeanResponseTimesFile # write_sim_files

# move up initialisation functions
export
	initCompTable!, initDmexclp!, initPriorityList!, initZhangIp!, initTemp0!, initTemp1!, initTemp2!

# misc functions
export
	printEvent, # event
	findNearestNode, # graph
	findNearestNodeInGrid, # grid
	isSameLocation, squareDist, normDist, offRoadTravelTime, linearInterpLocation, randLocation, # location
	isFNodeInRGraph, shortestPathNextRNode, shortestPathNextRArc, shortestPathTravelTime, shortestPath, findRArcFromFNodeToFNode, # network
	rasterRandLocations, printRasterSize, # raster
	printSimStats, printAmbsStats, printCallsStats, calcBatchMeans, calcBatchMeanResponseTimes, meanErrorPlot, calcAR0DurbinWatsonTestPValue, # statistics
	makeRandDeploymentPolicy, makeRandDeploymentPolicies, setAmbStation!, applyDeploymentPolicy!, simulateDeploymentPolicy!, simulateDeploymentPolicies! # deployment

# types
export
	Location, Node, Arc, Graph, NetTravel, Network, TravelMode, Travel, Route,
	Event, Ambulance, Call, Hospital, Station,
	Map, GridSearchRect, GridRect, Grid, Raster, RasterSampler, DemandMode, Demand,
	CompTableData, DmexclpData, PriorityListData, ZhangIpData, Temp0Data, Temp1Data, Temp2Data, MoveUpData,
	File, Table, Resimulation, Simulation,
	DistrRng

# defs - consts
export
	Float, FloatSpTime, IntRNode, IntFadj, # type alias
	nullIndex, nullX, nullY, nullTime, nullDist, # nulls
	Depol

# defs - enums
export
	Priority, nullPriority, highPriority, medPriority, lowPriority, # priorities
	AmbClass, nullAmbClass, als, bls, # ambulance classes
	EventForm, nullEvent, ambGoesToSleep, ambWakesUp, callArrives, considerDispatch, ambDispatched, ambReachesCall, ambGoesToHospital, ambReachesHospital, ambBecomesIdle, ambReachesStation, ambRedirected, considerMoveUp, ambMoveUp,
	AmbStatus, ambNullStatus, ambSleeping, ambIdleAtStation, ambGoingToCall, ambAtCall, ambGoingToHospital, ambAtHospital, ambGoingToStation,
	CallStatus, callNullStatus, callScreening, callQueued, callWaitingForAmb, callOnSceneCare, callGoingToHospital, callAtHospital, callProcessed,
	MoveUpModule, nullMoveUpModule, compTableModule, dmexclpModule, priorityListModule, zhangIpModule, temp0Module, temp1Module, temp2Module

include("defs.jl")

include("misc/rand.jl")

include("types/types.jl")
include("types/call.jl")
include("types/demand.jl")
include("types/event.jl")
include("types/graph.jl")
include("types/grid.jl")
include("types/location.jl")
include("types/network.jl")
include("types/raster.jl")
include("types/route.jl")
include("types/travel.jl")

include("file/file_io.jl")
include("file/gen_sim_files.jl")
include("file/read_sim_files.jl")
include("file/write_sim_files.jl")

include("decision/call.jl")
include("decision/dispatch.jl")

include("decision/move_up/comp_table.jl")
include("decision/move_up/dmexclp.jl")
include("decision/move_up/move_up_common.jl")
include("decision/move_up/priority_list.jl")
include("decision/move_up/temp0.jl")
include("decision/move_up/temp1.jl")
include("decision/move_up/temp2.jl")
include("decision/move_up/zhang_ip.jl")

include("resimulation.jl")
include("run_config.jl")
include("simulation.jl")
include("statistics.jl")

include("animation/animation.jl")

include("misc/deployment.jl")

end
