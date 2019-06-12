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

__precompile__()
module JEMSS

# animation
using HTTP
using Sockets
using WebSockets
using JSON

# files
using DelimitedFiles
using LightXML
using Serialization

# misc
using Base64
using CRC32c
using Mmap
using Pkg
using Printf
using Random
import Random: AbstractRNG, GLOBAL_RNG, MersenneTwister
using SparseArrays
using LightGraphs
using ArchGDAL
using Hungarian

# optimisation
using JuMP
using Cbc
using GLPK
using GLPKMathProgInterface
if haskey(Pkg.installed(), "Gurobi") && haskey(ENV, "GUROBI_HOME") using Gurobi end

# statistics
using Distributions
using HypothesisTests
using StatsBase
using StatsFuns

# simulation functions
export
	initSim, runConfig, # run_config
	simulate!, simulateToTime!, simulateToEnd!, backup!, reset!, simulateNextEvent!, # simulation
	animate!, animate, # animation
	backupSim!, resetSim! # compat

# file functions
export
	readDlmFileNextLine!, readDlmFile, openNewFile, writeDlmLine!, arrayDict, writeTablesToFile!, writeTablesToFile, readTablesFromFile, readTablesFromData, tableRowsFieldDicts, fileChecksum, serializeToFile, deserializeFile, joinPathIfNotAbs, interpolateString, xmlFileRoot, findElt, eltContent, eltContentVal, eltContentInterpVal, childrenNodeNames, selectXmlFile, # file_io
	runGenConfig, # gen_sim_files
	readAmbsFile, readArcsFile, readCallsFile, readCompTableFile, readDemandFile, readDemandCoverageFile, readEventsFile, readGeoFile, readHospitalsFile, readMapFile, readNodesFile, readPrioritiesFile, readPriorityListFile, readRasterFile, readRedispatchFile, readRNetTravelsFile, readStationsFile, readTravelFile, readDeploymentsFile, readZhangIpParamsFile, # read_sim_files
	writeAmbsFile, writeArcsFile, writeCallsFile, writeDemandFile, writeDemandCoverageFile, writeHospitalsFile, writeMapFile, writeNodesFile, writePrioritiesFile, writeRedispatchFile, writeRNetTravelsFile, writeSimPeriodStatsListFile, writeStationsFile, writeTravelFile, openOutputFiles!, closeOutputFiles!, writeEventToFile!, writeStatsFiles!, writeDeploymentsFile, writeBatchMeanResponseTimesFile # write_sim_files

# move up initialisation functions
export
	initCompTable!, initDmexclp!, initPriorityList!, initZhangIp!, initTemp0!, initTemp1!, initTemp2!

# misc functions
export
	initDemand!, initDemandCoverage!, # demand
	printEvent, # event
	findNearestNode, # graph
	findNearestNodeInGrid, # grid
	isSameLocation, squareDist, normDist, offRoadTravelTime, linearInterpLocation, randLocation, # location
	isFNodeInRGraph, shortestPathNextRNode, shortestPathNextRArc, shortestPathData, shortestPathTravelTime, shortestPath, findRArcFromFNodeToFNode, # network
	rasterRandLocations, printRasterSize, # raster
	shortestRouteTravelTime!, # route
	getCallResponseTimes, getAvgCallResponseTime, getCallsReachedInTime, countCallsReachedInTime, printSimStats, printAmbsStats, printCallsStats, printHospitalsStats, calcBatchMeans, calcBatchMeanResponseTimes, meanErrorPlot, calcAR0DurbinWatsonTestPValue, # statistics
	checkCompTable, checkCompTableIsNested, nestCompTable, unnestCompTable, makeRandNestedCompTable, # compliance table
	makeRandDeployment, makeRandDeployments, deploymentToStationsNumAmbs, stationsNumAmbsToDeployment, getDeployment, getStationsNumAmbs, setAmbStation!, applyDeployment!, applyStationsNumAmbs!, simulateDeployment!, simulateDeployments!, # deployment
	checkPriorityList, makeRandPriorityList, # priority list
	solveMexclp! # mexclp

# types
export
	Location, Point, Node, Arc, Graph, NetTravel, Network, TravelMode, Travel, Route,
	Event, Ambulance, Call, Hospital, Station, Redispatch,
	Map, GridSearchRect, GridRect, Grid, Raster, RasterSampler, DemandMode, Demand, PointsCoverageMode, DemandCoverage,
	CompTableData, DdsmData, DmexclpData, PriorityListData, ZhangIpData, Temp0Data, Temp1Data, Temp2Data, MoveUpData,
	AmbulanceStats, CallStats, HospitalStats, StationStats, SimPeriodStats, SimStats,
	File, Table, Resimulation, Simulation,
	DistrRng

# defs - consts
export
	Float, FloatSpTime, IntRNode, IntFadj, # type alias
	nullIndex, nullX, nullY, nullTime, nullDist, # nulls
	priorities, numPriorities, # priorities
	Deployment, CompTable, NestedCompTable, PriorityList

# defs - enums
export
	Priority, nullPriority, highPriority, medPriority, lowPriority, # priorities
	AmbClass, nullAmbClass, als, bls, # ambulance classes
	EventForm, nullEvent, ambGoesToSleep, ambWakesUp, callArrives, considerDispatch, ambDispatched, ambReachesCall, ambGoesToHospital, ambReachesHospital, ambBecomesIdle, ambReachesStation, ambRedirected, considerMoveUp, ambMoveUp,
	AmbStatus, ambNullStatus, ambSleeping, ambIdleAtStation, ambGoingToCall, ambAtCall, ambGoingToHospital, ambAtHospital, ambGoingToStation,
	CallStatus, callNullStatus, callScreening, callQueued, callWaitingForAmb, callOnSceneCare, callGoingToHospital, callAtHospital, callProcessed,
	MoveUpModule, nullMoveUpModule, compTableModule, dmexclpModule, priorityListModule, zhangIpModule, temp0Module, temp1Module, temp2Module

# deprecated
export
	Depol, makeRandDeploymentPolicy, makeRandDeploymentPolicies, applyDeploymentPolicy!, simulateDeploymentPolicy!, simulateDeploymentPolicies!, writeDeploymentPoliciesFile, readDeploymentPoliciesFile # renamed "deployment policy" to "deployment"

include("defs.jl")

include("misc/stream.jl")
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
include("types/statistics.jl")
include("types/travel.jl")

include("file/file_io.jl")
include("file/gen_sim_files.jl")
include("file/read_sim_files.jl")
include("file/write_sim_files.jl")

include("decision/call.jl")
include("decision/dispatch.jl")

include("decision/move_up/comp_table.jl")
include("decision/move_up/ddsm.jl")
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

include("optim/mexclp.jl")

include("deprecated.jl")

end
