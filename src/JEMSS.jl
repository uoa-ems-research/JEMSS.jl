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
using OffsetArrays
using Base.Iterators: Stateful
using Parameters
using Base.Threads

# optimisation
using JuMP
using Cbc
using GLPK
using Gurobi

# statistics
using Distributions
using HypothesisTests
using Statistics
using StatsBase
using StatsFuns

# simulation functions
export
    initSim, runConfig, # run_config
    simulate!, simulateToTime!, simulateToEnd!, backup!, reset!, simulateNextEvent!, # simulation
    setSimReps!, simulateRep!, simulateReps!, resetRep!, resetReps!, makeRepsRunnable!, # replication
    animate!, animate, # animation
    backupSim!, resetSim! # compat

# file functions
export
    readDlmFileNextLine!, readDlmFile, openNewFile, writeDlmLine!, arrayDict, writeTablesToFile!, writeTablesToFile, readTablesFromFile, readTablesFromData, tableRowsFieldDicts, fileChecksum, serializeToFile, deserializeFile, joinPathIfNotAbs, interpolateString, xmlFileRoot, findElt, eltContent, eltContentVal, eltContentInterpVal, childrenNodeNames, selectXmlFile, # file_io
    readGenConfig, runGenConfig, runGenConfigCalls, makeCalls, # gen_sim_files
    readAmbsFile, readArcsFile, readCallsFile, readCompTableFile, readDemandFile, readDemandCoverageFile, readEventsFile, readGeoFile, readHospitalsFile, readMapFile, readMobilisationDelayFile, readMultiCompTableFile, readNodesFile, readPrioritiesFile, readPriorityListFile, readPriorityListsFile, readRasterFile, readRedispatchFile, readRNetTravelsFile, readStationsFile, readStatsControlFile, readStatsDictFile, readTravelFile, readDeploymentsFile, readZhangIpParamsFile, # read_sim_files
    writeAmbsFile, writeArcsFile, writeCallsFile, writeDemandFile, writeDemandCoverageFile, writeHospitalsFile, writeMapFile, writeMultiCompTableFile, writeNodesFile, writePrioritiesFile, writePriorityListFile, writePriorityListsFile, writeRasterFile, writeRedispatchFile, writeRNetTravelsFile, writeStationsFile, writeTravelFile, writeZhangIpParamsFile, openOutputFiles!, writeOutputFiles, writeMiscOutputFiles, closeOutputFiles!, writeEventToFile!, writeDeploymentsFile, writeBatchMeanResponseDurationsFile, # write_sim_files
    writeStatsFiles, writeAmbsStatsFile, writeCallsStatsFile, writeHospitalsStatsFile, writeStationsStatsFile, writeStatsDictFile # write_sim_files - stats

# move up
export
    initCompTable!, initDmexclp!, initMultiCompTable!, initPriorityList!, initZhangIp!, initTemp0!, initTemp1!, initTemp2!,
    setMoveUpModule!

# misc functions
export
    isBusy, isFree, isWorking, isGoingToStation, isTravelling, # ambulance
    setSimCalls!, # call
    calcCoverBound!, initCoverBound, approxDistr, convolute, makeCdf, binarySearch, calcAmbBusyDurationProbUpperBounds, calcAmbBusyDurationLowerBoundDistrs!, calcNumAmbsMaxCoverageFrac, calcQueuedDurationsMaxCoverageFrac, initCoverBoundSim, simulateCoverBound!, simulateCoverBoundLowerBound!, writeAmbBusyDurationProbUpperBoundsFile, readAmbBusyDurationProbUpperBoundsFile, writeQueuedDurationsMaxCoverageFracFile, readQueuedDurationsMaxCoverageFracFile, # cover bound
    initDemand!, initDemandCoverage!, # demand
    printEvent, # event
    findNearestNode, # grid
    squareDist, normDist, offRoadTravelTime, linearInterpLocation, randLocation, # location
    isFNodeInRGraph, shortestPathNextRNode, shortestPathNextRArc, shortestPathData, shortestPathTravelTime, shortestPathDistance, shortestPath, findRArcFromFNodeToFNode, # network
    rasterRandLocations, printRasterSize, # raster
    shortestRouteTravelTime!, # route
    getCallResponseDurations, getAvgCallResponseDuration, getCallsReachedInTime, countCallsReachedInTime, printSimStats, printAmbsStats, printCallsStats, printHospitalsStats, calcBatchMeans, calcBatchMeanResponseDurations, meanErrorPlot, calcAR0DurbinWatsonTestPValue, tDistrHalfWidth, confInterval, getPeriodStatsList, getRepsPeriodStatsList, getRepPeriodStats, statsDictFromPeriodStatsList, # statistics
    checkCompTable, isCompTableNested, nestCompTable, unnestCompTable, makeRandNestedCompTable, # compliance table
    checkMultiCompTable, # multi-compliance table
    makeRandDeployment, makeRandDeployments, deploymentToStationsNumAmbs, stationsNumAmbsToDeployment, getDeployment, getStationsNumAmbs, setAmbStation!, applyDeployment!, applyStationsNumAmbs!, simulateDeployment!, simulateDeployments!, # deployment
    checkPriorityList, makeRandPriorityList, # priority list
    solveMclp!, solveMclp, # mclp
    solveMexclp!, # mexclp
    solveNestedMclp!, solveNestedMclp, # nested mclp
    solvePMedian, # p-median
    flatten, # dict
    runParallel! # parallel

# types
export
    Location, Point, Node, Arc, Graph, NetTravel, Network, TravelMode, Travel, Route,
    Event, Ambulance, Call, Hospital, Station, Redispatch,
    Map, GridSearchRect, GridRect, Grid, Raster, RasterSampler, DemandMode, Demand, PointsCoverageMode, DemandCoverage,
    CoverBoundSimRep, CoverBoundSim, CoverBoundMode, CoverBound,
    CompTableData, DdsmData, DmexclpData, MultiCompTableData, PriorityListData, ZhangIpData, Temp0Data, Temp1Data, Temp2Data, MoveUpData,
    MeanAndHalfWidth, AmbulanceStats, CallStats, HospitalStats, StationStats, SimPeriodStats, SimStats,
    File, Table, Resimulation, Simulation,
    DistrRng, GenConfig, MobilisationDelay

# defs - consts
export
    Float, FloatSpTime, FloatSpDist, IntRNode, IntFadj, # type alias
    nullIndex, nullX, nullY, nullTime, nullDist, nullHist, NullHist, # nulls
    priorities, numPriorities, # priorities
    Deployment, CompTable, MultiCompTable, NestedCompTable, PriorityList

# defs - enums
export
    Priority, nullPriority, highPriority, medPriority, lowPriority, # priorities
    AmbClass, nullAmbClass, als, bls, # ambulance classes
    EventForm, nullEvent, ambGoesToSleep, ambWakesUp, callArrives, considerDispatch, ambDispatched, ambMobilised, ambReachesCall, ambGoesToHospital, ambReachesHospital, ambBecomesFree, ambReturnsToStation, ambReachesStation, considerMoveUp, ambMoveUpToStation,
    AmbStatus, ambNullStatus, ambSleeping, ambIdleAtStation, ambMobilising, ambGoingToCall, ambAtCall, ambGoingToHospital, ambAtHospital, ambFreeAfterCall, ambReturningToStation, ambMovingUpToStation,
    AmbStatusSet, ambWorking, ambBusy, ambFree, ambTravelling, ambGoingToStation,
    CallStatus, callNullStatus, callScreening, callQueued, callWaitingForAmb, callOnSceneTreatment, callGoingToHospital, callAtHospital, callProcessed,
    RouteStatus, routeNullStatus, routeBeforeStartNode, routeOnPath, routeAfterEndNode,
    MoveUpModule, nullMoveUpModule, compTableModule, dmexclpModule, multiCompTableModule, priorityListModule, zhangIpModule, temp0Module, temp1Module, temp2Module

# deprecated
export
    Depol, makeRandDeploymentPolicy, makeRandDeploymentPolicies, applyDeploymentPolicy!, simulateDeploymentPolicy!, simulateDeploymentPolicies!, writeDeploymentPoliciesFile, readDeploymentPoliciesFile, # renamed "deployment policy" to "deployment"
    isSameLocation, findNearestNodeInGrid

include("defs.jl")

include("misc/dict.jl")
include("misc/histogram.jl")
include("misc/parallel.jl")
include("misc/rand.jl")
include("misc/stream.jl")

include("types/types.jl")
include("types/ambulance.jl")
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
include("types/station.jl")
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
include("decision/move_up/multi_comp_table.jl")
include("decision/move_up/priority_list.jl")
include("decision/move_up/temp0.jl")
include("decision/move_up/temp1.jl")
include("decision/move_up/temp2.jl")
include("decision/move_up/zhang_ip.jl")

include("resimulation.jl")
include("run_config.jl")
include("simulation.jl")
include("replication.jl")
include("statistics.jl")

include("animation/animation.jl")

include("misc/deployment.jl")

include("optim/mclp.jl")
include("optim/mexclp.jl")
include("optim/p_median.jl")
include("optim/cover_bound.jl")
include("optim/nested_mclp.jl")

include("deprecated.jl")

end
