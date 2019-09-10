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

# common definitions

const Float = Float64
# const Int = Int64
# const UInt = UInt64

# float and int types for storing all-pairs shortest paths data
# can reduce precision in order to reduce memory used
const FloatSpTime = Float64 # precision for storing shortest path travel times
const FloatSpDist = Float64 # precision for storing shortest path travel distances
const IntRNode = Int16 # precision for storing number of nodes in rGraph
const IntFadj = Int8 # precision for storing maximum number of nodes adjacent to any node in rGraph (= length of longest vector in network.rGraph.fadjList)

const sourceDir = @__DIR__
const sourcePath = sourceDir # for compat
const jemssDir = realpath(joinpath(sourceDir, ".."))

const pkgVersions = Pkg.installed()

const Deployment = Vector{Int} # deployment[i] gives index of station to deploy ambulance i to
const CompTable = Array{Int,2} # compTable[i,j] gives number of ambulances to place at station j when there are i free ambulances
const NestedCompTable = Vector{Int} # nested compliance table can be represented as a list of station indices; nestedCompTable[1:i] gives indices of stations that i free ambulances should be placed at
const PriorityList = Vector{Int} # a priority list gives the order of preference for which ambulances should be redeployed (moved up) to stations; priorityList[1] has the index of the station with highest priority

# run modes
const debugMode = false
const checkMode = true # for data checking, e.g. assertions that are checked frequently

# file chars
const delimiter = ','
const newline = "\r\n"

# misc null values
const nullIndex = -1
const nullX = -200 # outside lat/lon range
const nullY = -200
const nullTime = -1.0
const nullDist = -1.0
nullFunction() = nothing

# call/travel priorities
@enum Priority nullPriority=0 highPriority=1 medPriority=2 lowPriority=3
const priorities = setdiff([instances(Priority)...], [nullPriority])
const numPriorities = length(priorities)

# ambulance classes
@enum AmbClass nullAmbClass=0 als=1 bls=2 # als = advanced life support, bls = basic life support

# event types/forms
@enum EventForm nullEvent ambGoesToSleep ambWakesUp callArrives considerDispatch ambDispatched ambReachesCall ambGoesToHospital ambReachesHospital ambBecomesFree ambReturnsToStation ambReachesStation ambRedirected considerMoveUp ambMoveUpToStation

# ambulance statuses
@enum AmbStatus ambNullStatus ambSleeping ambIdleAtStation ambGoingToCall ambAtCall ambGoingToHospital ambAtHospital ambFreeAfterCall ambReturningToStation ambMovingUpToStation
@enum AmbStatusSet ambWorking ambBusy ambFree ambTravelling ambGoingToStation

# call statuses
@enum CallStatus callNullStatus callScreening callQueued callWaitingForAmb callOnSceneTreatment callGoingToHospital callAtHospital callProcessed

# route statuses
@enum RouteStatus routeNullStatus routeBeforeStartNode routeOnPath routeAfterEndNode

@enum MoveUpModule nullMoveUpModule compTableModule ddsmModule dmexclpModule priorityListModule zhangIpModule temp0Module temp1Module temp2Module
