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

mutable struct Location
	x::Float # longitude, or other
	y::Float # latitude, or other
	
	Location() = new(nullX, nullY)
	Location(x,y) = new(x,y)
end

mutable struct Point
	index::Int
	location::Location
	value::Any
	
	# node nearest to location
	nearestNodeIndex::Int
	nearestNodeDist::Float
	
	Point() = new(nullIndex, Location(), nothing,
		nullIndex, nullDist)
end

mutable struct Node
	index::Int
	location::Location
	offRoadAccess::Bool # if node can be used to get on-road and off-road
	
	fields::Dict{String,Any} # additional data, not used by simulation
	
	Node() = new(nullIndex, Location(), true,
		Dict{String,Any}())
end

mutable struct Arc
	index::Int
	fromNodeIndex::Int
	toNodeIndex::Int
	distance::Float # NaN if not set
	
	fields::Dict{String,Any} # additional data, not used by simulation
	
	Arc() = new(nullIndex, nullIndex, nullIndex, NaN,
		Dict{String,Any}())
end

# graph data for network (actually a digraph)
mutable struct Graph
	# parameters:
	isReduced::Bool
	nodes::Vector{Node}
	arcs::Vector{Arc}
	
	arcDists::Vector{Float} # arcDists[i] = arcs[i].distance
	
	light::LightGraphs.DiGraph
	fadjList::Vector{Vector{Int}} # shorthand for light.fadjlist; fadjList[i] gives nodes that are connected to node i by an outgoing arc from node i
	badjList::Vector{Vector{Int}} # shorthand for light.badjList; badjList[i] gives nodes that are connected to node i by an outgoing arc to node i
	
	# for full graph:
	nodePairArcIndex::SparseMatrixCSC{Int,Int} # for node indices i,j, nodePairArcIndex[i,j] should return arc index for node pair
	# cannot always use nodePairArcIndex for reduced graph, there may be more than one arc from node i to j, use spNodePairArcIndex instead
	
	Graph(isReduced::Bool) = new(isReduced, [], [],
		[],
		LightGraphs.DiGraph(), [], [],
		spzeros(Int, 0, 0))
end

# travel data for network
mutable struct NetTravel
	# parameters:
	isReduced::Bool
	modeIndex::Int
	arcTimes::Vector{Float} # arcTimes[i] gives travel time along arc i
	
	# for use in reduced graph:
	arcDists::Vector{Float} # arcDists[i] = distance of arc i; reference to Graph.arcDists, useful for checking that loaded rNetTravel is correct
	fadjList::Vector{Vector{Int}} # = rGraph.fadjList; useful for checking that loaded rNetTravel is correct
	spTimes::Array{FloatSpTime,2} # spTimes[i,j] = shortest path time between node i and j
	spDists::Array{FloatSpDist,2} # spDists[i,j] = shortest path distance between node i and j
	spFadjIndex::Array{IntFadj,2} # for shortest path from rNode i to j, spFadjIndex[i,j] gives the index (in fadjList[i], see Graph) of the successor rNode of i for this path
	spNodePairArcIndex::SparseMatrixCSC{Int,Int} # spNodePairArcIndex[i,j] = index of arc incident to nodes i,j, if it provides shortest travel time
	spFadjArcList::Vector{Vector{Int}} # spFadjArcList[i][j] gives index of rArc from rNode[i] to jth outgoing rNode of rNode[i] (i.e. rGraph.fadjList[i][j])
	
	# for use in full graph:
	fNodeToRNodeTime::Vector{Dict{Int,Float}} # fNodeToRNodeTime[i][j] gives time from fNode[i] to rNode[j], as long as rNode[j] is in fNodeToRNodes[i] (see type Network)
	fNodeFromRNodeTime::Vector{Dict{Int,Float}} # fNodeFromRNodeTime[i][j] gives time to fNode[i] from rNode[j], as long as rNode[j] is in fNodeFromRNodes[i] (see type Network)
	rArcFNodesTimes::Vector{Vector{Float}} # rArcFNodesTimes[i][k] gives travel time from rArc[i].fromNodeIndex to rArcFNodes[i][k]
	
	# useful object to/from node data, for full graph (see commonFNodes in Network):
	commonFNodeToFNodeTime::Array{Float,2} # commonFNodeToFNodeTime[i,j] gives shortest path time from commonFNodes[i] to fNode j
	fNodeToCommonFNodeTime::Array{Float,2} # fNodeToCommonFNodeTime[i,j] gives shortest path time from fNode i to commonFNodes[j]
	commonFNodeToFNodeDist::Array{Float,2} # commonFNodeToFNodeDist[i,j] gives shortest path distance from commonFNodes[i] to fNode j
	fNodeToCommonFNodeDist::Array{Float,2} # fNodeToCommonFNodeDist[i,j] gives shortest path distance from fNode i to commonFNodes[j]
	commonFNodeToFNodeRNodes::Array{Tuple{Int,Int},2} # commonFNodeToFNodeRNodes[i,j] gives start and end rNode for shortest path from commonFNodes[i] to fNode j
	fNodeToCommonFNodeRNodes::Array{Tuple{Int,Int},2} # fNodeToCommonFNodeRNodes[i,j] gives start and end rNode for shortest path from fNode i to commonFNodes[j]
	fNodeNearestHospitalIndex::Vector{Int} # fNodeNearestHospitalIndex[i] gives index of nearest hospital from fNode[i]
	
	NetTravel(isReduced::Bool) = new(isReduced, nullIndex, [],
		[], [], Array{FloatSpTime,2}(undef,0,0), Array{FloatSpDist,2}(undef,0,0), Array{IntFadj,2}(undef,0,0), spzeros(Int, 0, 0), [],
		[], [], [],
		Array{Float,2}(undef,0,0), Array{Float,2}(undef,0,0), Array{Float,2}(undef,0,0), Array{Float,2}(undef,0,0), Array{Tuple{Int,Int},2}(undef,0,0), Array{Tuple{Int,Int},2}(undef,0,0), [])
end

mutable struct Network
	fGraph::Graph # full graph
	rGraph::Graph # reduced graph
	
	# travel data for each graph:
	fNetTravels::Vector{NetTravel} # fNetTravels[i] gives fGraph travel data for travel mode i
	rNetTravels::Vector{NetTravel} # rNetTravels[i] gives rGraph travel data for travel mode i
	
	# for converting between full and reduced graph:
	rNodeFNode::Vector{Int} # rNodeFNode[i] returns index of rNode i in fGraph
	fNodeRNode::Vector{Int} # fNodeRNode[i] returns index of fNode i in rGraph, nullIndex if not in rGraph
	rArcFNodes::Vector{Vector{Int}} # rArcFNodes[i] gives array of fNode indices that belong to rArc i (ordered: fromNodeIndex -> toNodeIndex)
	fNodeRArcs::Vector{Vector{Int}} # fNodeRArcs[i] gives indices of rArcs that fNode[i] is on
	rArcFNodeIndex::Vector{Dict{Int,Int}} # rArcFNodeIndex[i][j] gives index that fNode[j] appears in rArc[i]; should be same as findall(rArcFNodes[i] .== j), so rArcFNodes[i][rArcFNodeIndex[i][j]] = j
	fNodeToRNodes::Vector{Vector{Int}} # fNodeToRNodes[i] gives indices of rNodes that fNode[i] can travel to, and is "near" ("near" meaning that the fNode is on an rArc incident to rNode, or fNode is rNode in rGraph; it is similar to adjacency)
	fNodeFromRNodes::Vector{Vector{Int}} # fNodeFromRNodes[i] gives indices of rNodes that can travel to fNode[i], and is "near"
	fNodeToRNodeNextFNode::Vector{Dict{Int,Int}} # fNodeToRNodeNextFNode[i][j] gives index of fNode after fNode[i] on path to rNode[j], where j is in fNodeToRNodes[i]. Needed to find path from an fNode to an rNode
	rArcFArcs::Vector{Vector{Int}} # rArcFArcs[i][j] gives the index of the jth fArc on rArcs[i]
	
	# distance
	fNodeToRNodeDist::Vector{Dict{Int,Float}} # fNodeToRNodeDist[i][j] gives distance from fNode[i] to rNode[j], as long as rNode[j] is in fNodeToRNodes[i]
	fNodeFromRNodeDist::Vector{Dict{Int,Float}} # fNodeFromRNodeDist[i][j] gives distance to fNode[i] from rNode[j], as long as rNode[j] is in fNodeFromRNodes[i]
	
	# fNodes that are common start/end points for travel (e.g, fNodes nearest to stations, hospitals):
	commonFNodes::Vector{Int} # list of common fNodes; = findall(isFNodeCommon)
	isFNodeCommon::Vector{Bool} # isFNodeCommon[i] = true if fNode i is common, false otherwise
	fNodeCommonFNodeIndex::Vector{Int} # fNodeCommonFNodeIndex[i] gives index of fNode i in commonFNodes, if isFNodeCommon[i] = true
	
	Network() = new(Graph(false), Graph(true),
		[], [],
		[], [], [], [], [], [], [], [], [],
		[], [],
		[], [], [])
end

# travel information, for on-road and off-road
mutable struct TravelMode
	index::Int
	
	offRoadSpeed::Float
	
	fNetTravel::NetTravel # reference to single fNetTravel in type Network
	rNetTravel::NetTravel # reference to single rNetTravel in type Network
	
	TravelMode() = new(nullIndex,
		nullTime,
		NetTravel(false), NetTravel(true))
end

# for storing travel modes, travel sets (a set of travel modes apply to each time period),
# and conditions for when to use each travel set/mode
mutable struct Travel
	numModes::Int # number of travel modes
	numSets::Int # number of sets of travel modes; each set may contain multiple travel modes e.g. for different combinations of travel priorities and ambulance classes. Different sets can overlap, by containing the same travel modes.
	
	modes::Vector{TravelMode}
	modeLookup::Array{Int,2} # modeLookup[i,j] = index of travel mode to use for travel set i, travel priority j. Change this variable according to modelling needs
	
	setsStartTimes::Vector{Float} # setsStartTimes[i] gives time at which travel set setsTimeOrder[i] should be started (setsStartTimes[i+1] or Inf gives end time)
	setsTimeOrder::Vector{Int} # setsTimeOrder[i] gives travel set to start using at time setsStartTimes[i]
	recentSetsStartTimesIndex::Int # index of most recently used value in setsStartTimes (and setsTimeOrder), should only ever increase in value
	
	Travel() = new(nullIndex, nullIndex,
		[], Array{Int,2}(undef,0,0),
		[], [], nullIndex)
end

# for storing ambulance routes
# routes include a path (on fGraph), and may start/end off the graph
mutable struct Route
	priority::Priority # travel priority
	travelModeIndex::Int
	
	# start and end locations and times
	startLoc::Location
	startTime::Float
	endLoc::Location
	endTime::Float
	
	# start and end nodes and times
	startFNode::Int # index of first fNode in route
	startFNodeTime::Float # time at which startFNode is reached
	startFNodeDist::Float # distance from startLoc to startFNode; can be different from result from normDist() if starting part way along an arc
	startRNode::Int # index of first rNode in route
	startRNodeTime::Float # time at which startRNode is reached
	endRNode::Int # index of last rNode in route
	endRNodeTime::Float # time at which endRNode is reached
	endFNode::Int # index of last fNode in route
	endFNodeTime::Float # time at which endFNode is reached
	
	firstRArc::Int # index of first rArc in route; = nullIndex if startFNode == endFNode
	
	## fields that vary throughout the route are below here
	## most of the route code assumes that these fields only change to move "forward" through the route
	
	status::RouteStatus
	recentUpdateTime::Float # time at which route was most recently updated
	
	# recent rArc visited ("recent" means it is from the most recent route update/query)
	recentRArc::Int # index of rArc recently visited
	recentRArcStartTime::Float # time that travel on rArc started / would have started if had travelled from first node of recentRArc (can start part way along the arc)
	recentRArcEndTime::Float # time that travel on rArc should end / would end if were travelling to last node of recentRArc (can end part way along the arc)
	
	# recent fNode visited on recentRArc
	recentRArcRecentFNode::Int # index of fNode recently visited on recentRArc, rArcFNodes[recentRArc][recentRArcRecentFNode] gives index of fNode in fGraph
	recentFNode::Int # = rArcFNodes[recentRArc][recentRArcRecentFNode]
	recentFNodeTime::Float # time that recentFNode was reached
	
	# next fNode to visit on recentRArc
	recentRArcNextFNode::Int # recentRArcRecentFNode + 1
	nextFNode::Int # = rArcFNodes[recentRArc][recentRArcNextFNode]
	nextFNodeTime::Float # time that nextFNode will be reached
	nextFNodeDist::Float # distance to reach nextFNode; this is only set after calling getRouteNextNodeDist!()
	
	Route() = new(nullPriority, nullIndex,
		Location(), nullTime, Location(), nullTime,
		nullIndex, nullTime, nullDist, nullIndex, nullTime, nullIndex, nullTime, nullIndex, nullTime,
		nullIndex,
		routeNullStatus, nullTime,
		nullIndex, nullTime, nullTime,
		nullIndex, nullIndex, nullTime,
		nullIndex, nullIndex, nullTime, nullDist)
end

mutable struct Event
	index::Int # index of event in list of events that have occurred (not for sim.eventList)
	parentIndex::Int # index of parent event
	form::EventForm # "type" is taken
	time::Float
	ambIndex::Int
	callIndex::Int
	stationIndex::Int
	
	Event() = new(nullIndex, nullIndex, nullEvent, nullTime, nullIndex, nullIndex, nullIndex)
end

mutable struct Ambulance
	index::Int
	status::AmbStatus
	stationIndex::Int
	callIndex::Int
	route::Route
	event::Event # next/current event, useful for deleting next event from eventList
	class::AmbClass # class/type of ambulance
	attributes::Dict{String,Any}
	# schedule::Schedule # to be added

	dispatchTime::Float # time at which amb was dispatched
	mobilisationTime::Float # time at which amb will mobilise if there is mobilisation delay; = nullTime if not mobilising
	
	# for animation:
	currentLoc::Location
	movedLoc::Bool
	destLoc::Location # destination
	
	# count statistics:
	numCallsTreated::Int # total number of calls that ambulance provided treatment on scene
	numCallsTransported::Int # total number of calls transported to hospital
	numDispatches::Int
	numDispatchesFromStation::Int # total number of dispatches while at station
	numDispatchesWhileMobilising::Int # total number of dispatches while mobilising
	numDispatchesOnRoad::Int # total number of dispatches while on road
	numDispatchesOnFree::Int # total number of dispatches after ambulance becomes free
	numRedispatches::Int # number of times that ambulance is redispatched from one call to another
	numMoveUps::Int
	numMoveUpsFromStation::Int
	numMoveUpsOnRoad::Int
	numMoveUpsOnFree::Int
	numMoveUpsReturnToPrevStation::Int # number of times that ambulance is told to return to the station from which it started move up (see moveUpFromStationIndex)
	
	# status statistics
	statusDurations::Dict{Union{AmbStatus,AmbStatusSet},Float} # duration spent in each status (or status set); only updated after setting status
	statusDistances::Dict{Union{AmbStatus,AmbStatusSet},Float} # distance travelled in each status (or status set); only updated after setting status
	statusTransitionCounts::Array{Int,2} # statusTransitionCounts[Int(status1), Int(status2)] gives number of times that amb transitioned from status1 to status2 (in AmbStatus enum)
	
	# for calculating statistics:
	statusSetTime::Float # time at which status was last set, even if set to same status value
	prevStatus::AmbStatus
	moveUpFromStationIndex::Int # index from which ambulance started move up; = nullIndex if not moving up from station
	
	# dispatch statistics:
	dispatchStartLocCounts::Dict{Location,Int} # dispatchStartLocCounts[loc] gives number of times that ambulance was dispatched from location 'loc'
	
	# # hospital and station statistics:
	# atHospitalDurations::Vector{Float} # atHospitalDurations[i] gives the duration that the ambulance spent at hospitals[i]
	# atHospitalCounts::Vector{Int} # atHospitalCounts[i] gives the number of times that the ambulance was at hospitals[i]
	# idleAtStationDurations::Vector{Int} # idleAtStationDurations[i] is the duration that the ambulance spent idle at stations[i]
	
	# # redispatch statistics
	# redispatchCounts::Array{Int,2} # redispatchCounts[Int(p1),Int(p2)] gives number of redispatches from call of priority p1 to call of priority p2
	
	Ambulance() = new(nullIndex, ambNullStatus, nullIndex, nullIndex, Route(), Event(), nullAmbClass, Dict(),
		nullTime, nullTime,
		Location(), false, Location(),
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		Dict(), Dict(), Array{Int,2}(undef,0,0),
		nullTime, ambNullStatus, nullIndex,
		Dict())
end

mutable struct Call
	index::Int
	priority::Priority
	transport::Bool # true if requires transport to hospital
	hospitalIndex::Int # hospital that call should be transported to. If hospitalIndex == nullIndex, will transport to nearest hospital. See also chosenHospitalIndex
	location::Location # where call occurs
	attributes::Dict{String,Any}
	
	# time/duration params:
	arrivalTime::Float # time at which call arrives
	dispatchDelay::Float # delay between call arrival and being ready to dispatch an ambulance
	onSceneDuration::Float # time spent at call location
	handoverDuration::Float # for hospital handover
	
	# node nearest to call location:
	nearestNodeIndex::Int
	nearestNodeDist::Float
	
	# call state
	status::CallStatus
	ambIndex::Int # index of responding ambulance
	
	# for animation:
	currentLoc::Location
	movedLoc::Bool
	
	# statistics:
	dispatchTime::Float # time at which final ambulance was dispatched (= arrivalTime + dispatchDelay, unless call queued/bumped)
	ambArrivalTime::Float # time at which ambulance arrives on-site
	hospitalArrivalTime::Float # time at which ambulance arrives at hospital
	numBumps::Int # total number of times that call gets bumped due to redispatch
	wasQueued::Bool # whether call was queued or not
	ambDispatchLoc::Location # location of responding ambulance at moment of dispatch
	ambStatusBeforeDispatch::AmbStatus # status of ambulance just before dispatch
	chosenHospitalIndex::Int # hospital (if any) that call was transported to
	
	# duration statistics:
	queuedDuration::Float # duration between being ready to dispatch (but no ambulances are free) and dispatching an ambulance
	bumpedDuration::Float # total duration spent waiting for ambulances that were later redirected away from call due to bumping
	waitingForAmbDuration::Float # total duration spent waiting for ambulances that have been dispatched, even if call bumped and ambulance did not reach call
	responseDuration::Float # duration between call arrival and ambulance arrival at call location
	ambGoingToCallDuration::Float # duration between dispatch of ambulance that responded and ambulance reaching call
	transportDuration::Float # duration for transporting call to hospital
	serviceDuration::Float # duration between call arrival and ambulance finishing with call
	
	# for calculating statistics:
	statusSetTime::Float # time at which status was last set, even if set to same status value
	
	Call() = new(nullIndex, nullPriority, true, nullIndex, Location(), Dict(),
		nullTime, nullTime, nullTime, nullTime,
		nullIndex, nullDist,
		callNullStatus, nullIndex,
		Location(), false,
		nullTime, nullTime, nullTime, 0, false, Location(), ambNullStatus, nullIndex,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		nullTime)
end

mutable struct Hospital
	index::Int
	location::Location
	attributes::Dict{String,Any}
	
	nearestNodeIndex::Int
	nearestNodeDist::Float
	
	# for statistics:
	numCalls::Int # total number of calls taken to hospital
	# numCallsTotalHandoverDuration::OffsetArray{Float,1,Vector{Float}} # numCallsTotalHandoverDuration[k] gives total duration that hospital had k calls for handover
	
	Hospital() = new(nullIndex, Location(), Dict(),
		nullIndex, nullDist,
		0)
end

mutable struct Station
	index::Int
	location::Location
	capacity::Int # maximum number of ambulances that station can hold
	attributes::Dict{String,Any}
	
	nearestNodeIndex::Int
	nearestNodeDist::Float
	
	# statistics:
	numIdleAmbsTotalDuration::OffsetArray{Float,1,Vector{Float}} # numIdleAmbsTotalDuration[k] gives total duration that station had k idle ambulances
	currentNumIdleAmbs::Int # current number of idle ambulances at station
	currentNumIdleAmbsSetTime::Float # time at currentNumIdleAmbs was last set (even if set to same value)
	
	# move up statistics:
	moveUpStartLocCounts::Dict{Location,Int} # locations (and counts) from which ambulances have started move up towards this station
	
	Station() = new(nullIndex, Location(), 0, Dict(),
		nullIndex, nullDist,
		OffsetVector(Float[],0), 0, nullTime,
		Dict())
end

mutable struct MobilisationDelay
	use::Bool # true if using mobilisation delay
	
	distrRng::DistrRng
	expectedDuration::Float
	
	MobilisationDelay() = new(false, DistrRng(Normal(0,0), seed = 0), 0.0)
	MobilisationDelay(use, distrRng, expectedDuration) = new(use, distrRng, expectedDuration)
end

# conditions that determine ambulance redispatch behaviour
mutable struct Redispatch
	allow::Bool # true if redispatch is allowed, false otherwise
	conditions::Array{Bool,2} # conditions[Int(p1),Int(p2)] returns true if ambulance may be redispatched from its current call (priority p1) to another call (priority p2). Change this var to match modelling needs.
	
	Redispatch() = new(false, fill(false,numPriorities,numPriorities))
	function Redispatch(::Type{Val{:default}})
		redispatch = Redispatch()
		redispatch.allow = true
		@assert(!any(redispatch.conditions)) # all values should be false so far
		redispatch.conditions[Int(lowPriority), Int(highPriority)] = true
		redispatch.conditions[Int(medPriority), Int(highPriority)] = true
		return redispatch
	end
end

mutable struct Map
	xMin::Float
	xMax::Float
	yMin::Float
	yMax::Float
	xRange::Float # xMax - xMin
	yRange::Float # yMax - yMin
	xScale::Float # convert delta(x) to distance
	yScale::Float # convert delta(y) to distance
	# xRangeDist::Float # xRange * xScale
	# yRangeDist::Float # yRange * yScale
	
	Map() = new(nullX, nullX, nullY, nullY, nullDist, nullDist, nullDist, nullDist)
	Map(xMin, xMax, yMin, yMax, xScale, yScale) = new(xMin, xMax, yMin, yMax, xMax - xMin, yMax - yMin, xScale, yScale)
end

# grid search rectangle, part of type Grid
# for keeping track of search progress,
# while looking for nearest node in grid
mutable struct GridSearchRect
	# distance of search rectangle borders from a given location
	xDist::Vector{Float}
	yDist::Vector{Float}
	
	# min and max range of grid indices to search next
	ixSearch::Vector{Int}
	iySearch::Vector{Int}
	
	# min and max grid indices of search rectangle, for each direction
	# contains indices of next search also
	ixSearched::Vector{Int}
	iySearched::Vector{Int}
	
	GridSearchRect() = new([nullDist, nullDist], [nullDist, nullDist],
		[nullIndex, nullIndex], [nullIndex, nullIndex],
		[nullIndex, nullIndex], [nullIndex, nullIndex])
	GridSearchRect(ix, iy) = new([nullDist, nullDist], [nullDist, nullDist],
		[ix, ix], [iy, iy],
		[ix, ix], [iy, iy])
end

# grid rectangle, part of type Grid
mutable struct GridRect
	nodeIndices::Vector{Int}
	
	GridRect() = new([])
end

# grid breaks map region into rectangles,
# each rectangle stores graph node indices (from full graph)
# this is for quickly finding the nearest node to a location
mutable struct Grid
	nx::Int # number of divisions in x direction
	ny::Int # number of divisions in y direction
	xRange::Float # size of divisions in x direction, = map.xRange / nx
	yRange::Float # size of divisions in y direction, = map.yRange / ny
	xRangeDist::Float # distance of divisions in x direction, = xRange * map.xScale
	yRangeDist::Float # distance of divisions in y direction, = yRange * map.yScale
	rects::Array{GridRect,2} # rectangles
	searchRect::GridSearchRect
	
	Grid() = new(nullIndex, nullIndex, nullDist, nullDist, nullDist, nullDist, Array{GridRect,2}(undef, 0, 0), GridSearchRect())
	function Grid(map::Map, nx, ny)
		grid = new(nx, ny, nullDist, nullDist, nullDist, nullDist, Array{GridRect,2}(undef, nx, ny))
		for i = 1:nx
			for j = 1:ny
				grid.rects[i,j] = GridRect()
			end
		end
		grid.xRange = map.xRange / nx
		grid.yRange = map.yRange / ny
		grid.xRangeDist = grid.xRange * map.xScale
		grid.yRangeDist = grid.yRange * map.yScale
		grid.searchRect = GridSearchRect()
		return grid
	end
end

# simple raster type
# stores data for rectangular grid of cells,
# cell (i,j) has centre x[i], y[j], and value z[i,j]
# requires at least 2 cells in each direction
mutable struct Raster
	# parameters:
	x::Vector{Float}
	y::Vector{Float}
	z::Array{Float,2} # z[i,j] corresponds with x[i], y[j]
	
	nx::Int # length(x)
	ny::Int # length(y)
	dx::Float # x step size
	dy::Float # y step size
	
	Raster() = new([], [], Array{Float,2}(undef,0,0),
		0, 0, 0.0, 0.0)
	function Raster(x, y, z)
		nx = length(x)
		ny = length(y)
		@assert(nx > 1)
		@assert(ny > 1)
		@assert((nx, ny) == size(z))
		dx = (maximum(x) - minimum(x))/(nx-1)
		dy = (maximum(y) - minimum(y))/(ny-1)
		@assert(dx > 0)
		@assert(dy > 0)
		return new(x, y, z, nx, ny, dx, dy)
	end
end

# for generating random locations within a raster,
# where raster z values are proportional to probability of location being
# within corresponding cell (i.e. z values are for categorical distribution),
# and locations are generated uniformly within the cell
mutable struct RasterSampler
	raster::Raster
	cellDistrRng::DistrRng # to generate index of raster cell
	cellLocRng::AbstractRNG # rng used when generating location within raster cell
	
	function RasterSampler(raster::Raster, cellRng::AbstractRNG, cellLocRng::AbstractRNG)
		@assert(all(raster.z .>= 0)) # otherwise probability values will be negative
		p = raster.z[:] / sum(raster.z) # pdf for z
		cellSampler = sampler(Categorical(p))
		cellDistrRng = DistrRng(cellSampler, cellRng)
		return new(raster, cellDistrRng, cellLocRng)
	end
	function RasterSampler(raster::Raster, cellSeed::Int, cellLocSeed::Int)
		cellRng = (cellSeed >= 0 ? MersenneTwister(cellSeed) : MersenneTwister(rand(UInt32)))
		cellLocRng = (cellLocSeed >= 0 ? MersenneTwister(cellLocSeed) : MersenneTwister(rand(UInt32)))
		return RasterSampler(raster, cellRng, cellLocRng)
	end
end

# to model call demand for a single priority
# can use multiple demand modes to model change in demand with time, using Demand type
mutable struct DemandMode
	# parameters:
	index::Int # mode index
	rasterIndex::Int # index of demand raster in Demand.rasters
	priority::Priority # demand priority
	arrivalRate::Float # demand per day
	
	raster::Raster # reference to Demand.rasters[rasterIndex]
	rasterMultiplier::Float # = arrivalRate / sum(raster.z), to scale raster.z values to match arrivalRate
	
	DemandMode() = new(nullIndex, nullIndex, nullPriority, nullTime,
		Raster(), 0.0)
end

# To model demand (calls).
# Stores demand modes, demand sets (a set of demand modes apply to each time period),
# and conditions for when to use each demand set/mode.
mutable struct Demand
	initialised::Bool
	numRasters::Int # number of demand rasters
	numModes::Int # number of demand modes
	numSets::Int # number of sets of demand modes; each set may contain multiple demand modes e.g. for different combinations of call priorities and ambulance classes. Different sets can overlap, by containing the same demand modes.
	
	rasters::Vector{Raster} # demand rasters, to model spatial demand. Raster cell z values are not actual demand, but are proportional to demand
	rasterFilenames::Vector{String} # rasterFilenames[i] gives filename of rasters[i]
	
	modes::Vector{DemandMode}
	modeLookup::Array{Int,2} # modeLookup[i,j] = index of demand mode to use for demand set i, demand priority j. Change this variable according to modelling needs
	
	setsStartTimes::Vector{Float} # setsStartTimes[i] gives time at which demand set setsTimeOrder[i] should be started (setsStartTimes[i+1] or Inf gives end time)
	setsTimeOrder::Vector{Int} # setsTimeOrder[i] gives demand set to start using at time setsStartTimes[i]
	recentSetsStartTimesIndex::Int # index of most recently used value in setsStartTimes (and setsTimeOrder), should only ever increase in value
	
	Demand() = new(false, nullIndex, nullIndex, nullIndex,
		[], [],
		[], Array{Int,2}(undef,0,0),
		[], [], nullIndex)
end

# For a given set of points, coverage time, and travel mode,
# for each station store the points that can be reached by
# travelling from the station within the coverage time.
mutable struct PointsCoverageMode
	# parameters:
	index::Int
	points::Vector{Point} # reference to DemandCoverage.points
	coverTime::Float # any demand that can be reached within this time gets additional coverage of 1, otherwise 0. Usually equal to target cover time minus dispatch delay.
	travelMode::TravelMode # reference to a travel mode
	
	pointSets::Vector{Vector{Int}} # pointSets[i] has set of all point indices covered by the same unique set of stations
	stationSets::Vector{Vector{Int}} # stationSets[i] = unique set of stations for pointSets[i]
	stationsCoverPointSets::Vector{Vector{Int}} # stationsCoverPointSets[i] = indices of pointSets covered by station i
	
	PointsCoverageMode() = new(nullIndex, [], nullTime, TravelMode(),
		[], [], [])
end

mutable struct DemandCoverage
	# params
	coverTimes::Dict{Priority,Float} # coverTimes[p] gives the target cover time for demand of priority p
	rasterCellNumRows::Int # number of rows of points to create per demand raster cell
	rasterCellNumCols::Int # number of columns of points to create per demand raster cell
	
	initialised::Bool
	points::Vector{Point} # demand is aggregated to points, same points are used for all demand rasters
	nodesPoints::Vector{Vector{Int}} # nodesPoints[i] gives indices of points for which node i is the nearest node
	rastersPointDemands::Vector{Vector{Float}} # rastersPointDemands[i][j] is demand at points[j] for Demand.rasters[i]
	
	pointsCoverageModes::Vector{PointsCoverageMode}
	pointsCoverageModeLookup::Vector{Dict{Float,Int}} # pointsCoverageModeLookup[TravelMode.index][coverTime] gives index of PointsCoverageMode
	pointSetsDemands::Array{Vector{Float},2} # pointSetsDemands[PointsCoverageMode.index, DemandMode.rasterIndex] gives relative demand values for each point set in PointsCoverageMode.pointSets, for Demand.rasters[rasterIndex]. Note that this needs to be multiplied by DemandMode.rasterMultiplier to get absolute (instead of relative) demand values.
	
	DemandCoverage() = new(Dict(), 0, 0,
		false, [], [], [],
		[], [], Array{Vector{Float},2}(undef,0,0))
	
	function DemandCoverage(coverTimes::Dict{Priority,Float}, rasterCellNumRows::Int, rasterCellNumCols::Int)
		dc = demandCoverage = DemandCoverage()
		dc.coverTimes = coverTimes
		dc.rasterCellNumRows = rasterCellNumRows
		dc.rasterCellNumCols = rasterCellNumCols
		return demandCoverage
	end
	
	function DemandCoverage(demandCoverage::DemandCoverage)
		dc = demandCoverage # shorthand
		return DemandCoverage(dc.coverTimes, dc.rasterCellNumRows, dc.rasterCellNumCols)
	end
end

# Mean and half the width of the confidence interval of the mean.
# Half-width is also known as half-length, or margin of error.
struct MeanAndHalfWidth
	mean::Float
	halfWidth::Float
	MeanAndHalfWidth(m, hw) = new(convert(Float,m), convert(Float,hw))
end

# cover bound simulation replication
mutable struct CoverBoundSimRep
	numCalls::Int
	numFreeAmbsCount::Vector{Int} # numFreeAmbsCount[i] = number of times that calls arrived while there were i free ambs
	queuedCallDurations::Vector{Float} # for calls that arrived when there were no free ambs, record how long each call was queued
	bound::Float # upper bound on fraction of in time responses
	
	CoverBoundSimRep() = new(0, [], [], 0)
end

# cover bound simulation
mutable struct CoverBoundSim
	numReps::Int
	numAmbs::Int
	
	warmUpDuration::Float
	minLastCallArrivalTime::Float
	
	reps::Vector{CoverBoundSimRep}
	bound::MeanAndHalfWidth # upper bound on fraction of in time responses
	
	interarrivalTimeDistrRng::DistrRng
	ambBusyDurationRng::AbstractRNG
	
	CoverBoundSim() = new(0, 0,
		nullTime, nullTime,
		[], MeanAndHalfWidth(0,0))
end

# cover bound data for a single combination of conditions (e.g. single demand mode and single travel mode)
mutable struct CoverBoundMode
	# parameters
	index::Int
	responseTravelMode::TravelMode # travel mode for ambulance response
	# demandPriorities::Vector{Priority} # all demand priorities that use this cover bound mode
	
	# travel times
	stationsToNodesTimes::Array{Float,2} # stationsToNodesTimes[i,j] = travel time from station i to node j
	stationsToPointsTimes::Array{Float,2} # stationsToPointsTimes[i,j] = travel time from station i to demand point j
	pointsToHospitalTimes::Vector{Float} # pointsToHospitalTimes[i] = travel time from point i to nearest hospital
	
	# group demand into sets
	# use in conjunction with demandCoverage, e.g. demandCoverage.nodesPoints
	nodeSets::Vector{Vector{Int}} # nodeSets[i] gives node indices that have same order in which ambulances can reach from each station
	nodeSetsStationList::Vector{Vector{Int}} # nodeSetsStationList[i] has list of stations ordered (from closest in time to furthest) from nodeSets[i]
	
	# distributions related to service
	distrs::Dict{String,Distribution} # distrs[name] gives distribution; name from Set("onSceneDuration", "handoverDuration")
	approxDistrs::Dict{String,Distribution} # approxDistrs[name] gives lower bound approximation (discrete, non-parametric distribution) of distrs[name]
	# cdfs::Dict{String,Vector{Float}} # to add? for "d1" and "d2"
	transportProb::Float # probability that call requires transport to hospital; assume call is taken to nearest hospital
	
	CoverBoundMode() = new(nullIndex, TravelMode(),
		Array{Float,2}(undef,0,0), Array{Float,2}(undef,0,0), [],
		[], [],
		Dict(), Dict(), 1.0)
end

# cover bound - a bound on optimal move-up performance
mutable struct CoverBound
	modes::Vector{CoverBoundMode}
	modeLookup::Vector{Int} # modeLookup[i] gives index of mode (in modes) to use for travel mode i. Will assume for now that there is only one demand set (demand is static) and one travel set (i.e. travel times are static).
	
	# parameters
	ambBusyDurationsToSample::Vector{Float} # will find upper bound on probability of amb being busy for each duration in this list
	accountForQueuedDurations::Bool # for queued calls, calculate maximum coverage according to time remaining to reach call on time; set to false for original cover bound which calculates coverage of queued calls as if the call was queued for no time
	queuedDurationsToSample::Vector{Float} # will solve MCLP for different coverage times based on how long the calls were queued; use if accountForQueuedDurations == true
	
	dispatchDelay::Float # assumed constant for simpler modelling; add any mobilisation delay to this
	
	# amb busy duration
	ambBusyDurationProbUpperBounds::Array{Float,2} # ambBusyDurationProbUpperBounds[i,j] = upper bound on probability of dispatched amb being busy for duration <= ambBusyDurationsToSample[i] for j free ambulances
	ambBusyDurationLowerBoundDistrs::Vector{Sampleable} # ambBusyDurationLowerBoundDistrs[i] is sampleable distribution of amb busy duration for i free ambs
	
	# coverage
	numAmbsMaxCoverageFrac::Vector{Float} # numAmbsMaxCoverageFrac[i] gives maximum fraction of demand that can be covered with i free ambulances (from 1:numStations), from solving MCLP
	queuedDurationsMaxCoverageFrac::Vector{Float} # queuedDurationsMaxCoverageFrac[i] = maximum coverage for call that was queued for a duration of queuedDurationsToSample[i]; use if accountForQueuedDurations == true
	
	sim::CoverBoundSim
	
	CoverBound() = new([], [],
		[], false, [],
		0.0,
		Array{Float,2}(undef,0,0), [],
		[], [],
		CoverBoundSim())
end

# move up data types
abstract type MoveUpDataType end
mutable struct EmptyMoveUpData <: MoveUpDataType end

# compliance table data
mutable struct CompTableData <: MoveUpDataType
	# parameters:
	compTable::CompTable # compTable[i,j] = number of ambulances to place at station j, with i total free ambs
	
	compTableStationSlots::Vector{Vector{Int}} # sum(compTableStationSlots[i] .== j) == compTable[i,j]
	
	# arrays for recycling:
	ambMovable::Vector{Bool} # ambMovable[i] = true if ambulance i can be moved up, false otherwise
	
	CompTableData() = new(CompTable(undef,0,0),
		[],
		[])
end

# ddsm - dynamic double standard model
mutable struct DdsmData <: MoveUpDataType
	# parameters:
	coverFractionTargetT1::Float # fraction of demand to cover within coverTimes[1]
	travelTimeCost::Float
	slackWeight::Float
	coverTimeDemandPriorities::Vector{Priority} # coverTimeDemandPriorities[i] is demand priority for coverTimes[i]
	
	coverTimes::Vector{Float} # target coverage response durations; coverTimes[1] <= coverTimes[2]
	options::Dict{Symbol,Any}
	
	DdsmData() = new(0.0, 0.0, 0.0, [],
		[], Dict())
end

# dmexclp - dynamic maximum expected coverage location problem
mutable struct DmexclpData <: MoveUpDataType
	# parameters:
	busyFraction::Float # fraction for which each ambulance is busy, approximate
	demandWeights::Dict{Priority,Float} # weight of each demand priority on the objective function
	# some other relevant parameters are stored in sim: demand, demandCoverage, responseTravelPriorities
	
	marginalBenefit::Vector{Float} # marginalBenefit[i] = benefit of adding an ith ambulance to cover single demand, calculated from busyFraction
	
	# arrays for recycling:
	stationNumFreeAmbs::Vector{Int} # stationNumFreeAmbs[i] = number of free ambulances assigned to station i
	stationMarginalCoverages::Vector{Float} # stationMarginalCoverages[i] gives extra coverage provided from placing newly freed ambulance at station i
	# pointSetsCoverCounts::Vector{Vector{Int}} # pointSetsCoverCounts[i][j] = number of free ambulances covering node set j, for demand.pointsCoverageModes i
	
	DmexclpData() = new(0.0, Dict(),
		[],
		[], [])
end

# multi-compliance table data
mutable struct MultiCompTableData <: MoveUpDataType
	# parameters:
	multiCompTable::MultiCompTable # multiCompTable[i][j][k] is jth option for number of ambs to assign to station k with i free ambs
	
	MultiCompTableData() = new(MultiCompTable())
end

mutable struct PriorityListData <: MoveUpDataType
	# parameters:
	priorityList::Vector{Int} # priorityList[i] gives station index that the ith free ambulance should be moved to
	
	# arrays for recycling:
	stationNumFreeAmbs::Vector{Int} # number of free ambulances assigned to each station
	
	PriorityListData() = new([],
		[])
end

mutable struct ZhangIpData <: MoveUpDataType
	# parameters:
	marginalBenefits::Vector{Vector{Float}} # marginalBenefits[i][j] = benefit of adding a jth ambulance to station i
	stationCapacities::Vector{Int} # stationCapacities[i] is the number of ambulances that station i can hold on completion of move up, should be <= stations[i].capacity
	travelTimeCost::Float # travel time cost multiplier
	onRoadMoveUpDiscountFactor::Float # discount of travel cost of move up of ambulance on-road and with "regret" travel time <= regretTravelTimeThreshold
	regretTravelTimeThreshold::Float
	expectedHospitalHandoverDuration::Float
	
	stationSlots::Vector{Int}
	benefitSlots::Vector{Float}
	marginalBenefitsDecreasing::Bool # true if benefit of adding ambulance j to a station is < benefit of adding ambulance j-1, for all stations
	stationSlotsOrderPairs::Array{Int,2} # stationSlotsOrderPairs[i,1:2] gives two stationSlots indices, first should be filled (with ambulance) before second
	
	ZhangIpData() = new([], [], 1.0, 1.0, nullTime, nullTime,
		[], [], false, Array{Int,2}(undef,0,0))
end

mutable struct Temp0Data <: MoveUpDataType
	# parameters:
	busyFraction::Float # ambulance busy fraction - should remove this, make marginalBenefit a parameter
	travelTimeCost::Float # travel time cost multiplier
	maxIdleAmbTravelTime::Float # max travel time for idle ambulances. 0.021 (days) is about 30 minutes
	maxNumNearestStations::Int # number of nearest stations to consider for each ambulance (may include ambulance's station)
	
	marginalBenefit::Vector{Float} # marginalBenefit[i] = benefit of adding an ith ambulance to a station
	
	Temp0Data() = new(0.0, 0.0, Inf, 0,
		[])
end

mutable struct Temp1Data <: MoveUpDataType
	# parameters:
	benefit::Array{Vector{Float},2} # benefit[i1,i2][k] = benefit from having k ambulances collectively at stations i1 and i2
	busyFraction::Float # ambulance busy fraction - should remove this, make marginalBenefit a parameter
	stationPairs::Vector{Vector{Int}} # stationPairs[i] gives ith station pair, which is a pair of station indices
	travelTimeCost::Float # travel time cost multiplier
	maxIdleAmbTravelTime::Float # max travel time for idle ambulances
	maxNumNearestStations::Int # number of nearest stations to consider for each ambulance (may include ambulance's station)
	
	marginalBenefit::Array{Vector{Float},2} # marginal benefit values from ambulance to station allocations, calculated from 'benefit'
	
	Temp1Data() = new(Array{Vector{Float},2}(undef,0,0), 0.0, Vector{Vector{Int}}(), 0.0, Inf, 0,
		Array{Vector{Float},2}(undef,0,0))
end

mutable struct Temp2Data <: MoveUpDataType
	# parameters:
	benefit::Array{Array{Float,2},2} # benefit[i1,i2][k1,k2] = benefit from having k1 and k2 ambulances at stations i1 and i2, respectively
	busyFraction::Float # ambulance busy fraction - should remove this, make marginalBenefit a parameter
	stationPairs::Vector{Vector{Int}} # stationPairs[i] gives ith station pair, which is a pair of station indices
	travelTimeCost::Float # travel time cost multiplier
	maxIdleAmbTravelTime::Float # max travel time for idle ambulances
	maxNumNearestStations::Int # number of nearest stations to consider for each ambulance (may include ambulance's station)
	
	marginalBenefit::Array{Array{Float,2},2} # marginal benefit values from ambulance to station allocations, calculated from 'benefit'
	
	Temp2Data() = new(Array{Array{Float,2},2}(undef,0,0), 0.0, Vector{Vector{Int}}(), 0.0, Inf, 0,
		Array{Array{Float,2},2}(undef,0,0))
end

mutable struct MoveUpData
	# parameters:
	useMoveUp::Bool
	moveUpModule::MoveUpModule # indicates move up module to be used
	
	# move up data types:
	compTableData::CompTableData
	ddsmData::DdsmData
	dmexclpData::DmexclpData
	multiCompTableData::MultiCompTableData
	priorityListData::PriorityListData
	zhangIpData::ZhangIpData
	temp0Data::Temp0Data
	temp1Data::Temp1Data
	temp2Data::Temp2Data
	
	MoveUpData() = new(false, nullMoveUpModule,
		CompTableData(), DdsmData(), DmexclpData(), MultiCompTableData(), PriorityListData(), ZhangIpData(), Temp0Data(), Temp1Data(), Temp2Data())
end

mutable struct File
	name::String
	path::String # includes name
	iostream::IOStream
	checksum::UInt
	
	# File() = new()
	File() = new("", "", IOStream(""), 0)
end

mutable struct EventsFile
	io::IOStream # write simulation trace of events to this file
	eventFilter::Dict{EventForm,Bool} # eventFilter[eventForm] = true if saving events of eventForm to file, false otherwise
	
	EventsFile() = new(IOStream(""), Dict([e => true for e in instances(EventForm)]))
end

mutable struct Resimulation
	# parameters:
	use::Bool # true if resimulating (will follow event trace), false otherwise
	timeTolerance::Float
	
	events::Vector{Event} # note that events[i].index may not necessarily equal i (events[i].index >= i)
	eventsChildren::Vector{Vector{Event}} # eventsChildren[i] gives events that are children of sim event i; length(eventsChildren) may be less than number of events if events file does not include all events
	prevEventIndex::Int # index of previous event in events field
	
	eventFilter::Dict{EventForm,Bool} # eventFilter[eventForm] = true if resimulating with events of eventForm
	doDispatch::Bool # true if can resimulate dispatch decisions, i.e. eventFilter[ambDispatched] == true
	doMoveUp::Bool # true if can resimulate move-up decisions, i.e. eventFilter[ambMoveUpToStation] == true
	
	Resimulation() = new(false, 0.0,
		[], [], nullIndex,
		Dict(), false, false)
end

# statistics for a single ambulance, or multiple ambulances
mutable struct AmbulanceStats
	ambIndex::Int # for single ambulance
	
	numCallsTreated::Int # total number of calls that ambulance provided treatment on scene
	numCallsTransported::Int # total number of calls transported to hospital
	
	numDispatches::Int # all dispatches
	numDispatchesFromStation::Int # total number of dispatches while at station
	numDispatchesWhileMobilising::Int # total number of dispatches while mobilising
	numDispatchesOnRoad::Int # total number of dispatches while on road
	numDispatchesOnFree::Int # total number of dispatches after ambulance becomes free
	numRedispatches::Int # number of times that ambulance is redispatched from one call to another
	
	numMoveUps::Int
	numMoveUpsFromStation::Int
	numMoveUpsOnRoad::Int
	numMoveUpsOnFree::Int
	numMoveUpsReturnToPrevStation::Int
	
	statusDurations::Dict{Union{AmbStatus,AmbStatusSet},Float} # duration spent in each status (or status set)
	statusDistances::Dict{Union{AmbStatus,AmbStatusSet},Float} # distance travelled in each status (or status set)
	statusTransitionCounts::Array{Int,2} # statusTransitionCounts[Int(status1), Int(status2)] gives number of times that amb transitioned from status1 to status2
	
	AmbulanceStats() = new(nullIndex,
		0, 0,
		0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0,
		Dict(), Dict(), Array{Int,2}(undef,0,0))
end

# statistics for a single call, or multiple calls
mutable struct CallStats
	callIndex::Int # for single call
	numCalls::Int # for multiple calls
	
	numQueued::Int
	numBumped::Int # number of calls bumped once or more
	numBumps::Int # number of bumps from all calls
	numTransports::Int
	numResponsesInTime::Int # number of in time ambulance responses
	
	# durations
	totalDispatchDelay::Float # delay between call arrival and being ready to dispatch an ambulance
	totalOnSceneDuration::Float # duration spent at call location
	totalHandoverDuration::Float # for hospital handover
	totalQueuedDuration::Float # duration between being ready to dispatch (but no ambulances are free) and dispatching an ambulance
	totalBumpedDuration::Float # total duration spent waiting for ambulances that were later redirected away from call due to bumping
	totalWaitingForAmbDuration::Float # total duration spent waiting for ambulances that have been dispatched, even if call bumped and ambulance did not reach call
	totalResponseDuration::Float # duration between call arrival and ambulance arrival at call location
	totalAmbGoingToCallDuration::Float # duration of ambulance that responded to reach call
	totalTransportDuration::Float # duration for transporting call to hospital
	totalServiceDuration::Float # duration between call arrival and ambulance finishing with call
	
	# totalQueueLengthDuration::OffsetArray{Float,1,Vector{Float}} # totalQueueLengthDuration[k] = total duration spent with k queued calls
	responseDurationHist::Histogram
	
	CallStats() = new(nullIndex, 0,
		0, 0, 0, 0, 0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		NullHist())
end

# statistics for a single hospital, or multiple hospitals
mutable struct HospitalStats
	hospitalIndex::Int # for single hospital
	
	numCalls::Int # total number of calls taken to hospital
	
	HospitalStats() = new(nullIndex,
		0)
end

# statistics for single station, or multiple stations
mutable struct StationStats
	stationIndex::Int # for single station
	
	numIdleAmbsTotalDuration::OffsetArray{Float,1,Vector{Float}} # numIdleAmbsTotalDuration[k] gives total duration that station had k idle ambulances
	
	# numMoveUpsOnFree::Int # number of move-ups to this station of ambulances that have just become free
	# numMoveUpsOnRoad::Int # number of move-ups to this station that were from ambulances that were driving on the road
	# numMoveUpsFromHospitals::Vector{Int} # numMoveUpsFromHospitals[i] is the number of move-ups to this station of ambulances from hospitals[i]
	# numMoveUpsFromStations::Vector{Int} # numMoveUpsFromStations[i] is the number of move-ups of ambulances from stations[i] to this station
	# numMoveUpsToStations::Vector{Int} # numMoveUpsToStations[i] is the number of move-ups of ambulances from this station to stations[i]
	# also count how many move-ups were actually completed, and how many were cut short?
	
	StationStats() = new(nullIndex,
		OffsetVector(Float[],0))
		# 0, 0, [], [], [])
end

# simulation statistics for a period
mutable struct SimPeriodStats
	# statistics are for [startTime, endTime)
	startTime::Float # start time for stats
	endTime::Float # end time for stats
	duration::Float # = endTime - startTime
	
	ambulances::Vector{AmbulanceStats}
	hospitals::Vector{HospitalStats}
	stations::Vector{StationStats}
	
	# aggregate values
	ambulance::AmbulanceStats # for all ambulances
	call::CallStats # statistics of calls that arrived within [startTime, endTime)
	callPriorities::Dict{Priority,CallStats} # callPriorities[p] is statistics of priority p calls that arrived within [startTime, endTime)
	hospital::HospitalStats # for all hospitals
	station::StationStats # for all stations
	
	# numAmbsFreeTotalDuration::OffsetArray{Float,1,Vector{Float}} # numAmbsFreeTotalDuration[k] gives total duration that k ambulances were free (working but not busy)
	
	SimPeriodStats() = new(nullTime, nullTime, nullTime,
		[], [], [],
		AmbulanceStats(), CallStats(), Dict(), HospitalStats(), StationStats())
end

mutable struct SimStats
	captures::Vector{SimPeriodStats} # from sim start to some time
	periods::Vector{SimPeriodStats} # arbitrary start and end times; calculated from captures, periods[i] = captures[i] - captures[i-1]
	
	# params
	doCapture::Bool # true if capturing statistics during simulation
	warmUpDuration::Float # first capture will be at sim.startTime + warmUpDuration (unless warmUpDuration = 0)
	periodDurationsIter::Stateful # iterator that returns duration of next stats period
	nextCaptureTime::Float # time at which statistics will next be captured; updated according to periodDurationsIter
	
	# sim timestamps
	simStartTime::Float
	simEndTime::Float
	warmUpEndTime::Float # = simStartTime + warmUpDuration
	lastCallArrivalTime::Float # useful for ignoring cool-down period
	
	# misc
	recordDispatchStartLocCounts::Bool # if true, will save starting location of ambulance dispatches in ambulance.dispatchStartLocCounts
	recordMoveUpStartLocCounts::Bool # if true, will save starting location of ambulance move up in station.moveUpStartLocCounts
	
	# response duration histogram
	recordResponseDurationHist::Bool # if true, will save call response durations to histogram stored in callStats.responseDurationHist
	responseDurationHistBinWidth::Float # width (duration) of bins for callStats.responseDurationHist
	responseDurationHistNumBins::Int # number of bins for callStats.responseDurationHist
	
	SimStats() = new([], [],
		true, 0.0, Stateful([]), Inf,
		nullTime, nullTime, nullTime, nullTime,
		false, false,
		false, nullTime, 0)
end

mutable struct Simulation
	startTime::Float
	time::Float # time of most recent event, or time of recent animation frame if animating
	endTime::Float # calculated after simulating
	
	# world:
	net::Network
	travel::Travel
	map::Map
	grid::Grid
	
	ambulances::Vector{Ambulance}
	calls::Vector{Call}
	hospitals::Vector{Hospital}
	stations::Vector{Station}
	
	# shorthand:
	numAmbs::Int # length(ambulances)
	numCalls::Int # length(calls)
	numHospitals::Int # length(hospitals)
	numStations::Int # length(stations)
	
	eventList::Vector{Event} # events to occur now or in future
	eventIndex::Int # index of event in events that have occurred
	queuedCallList::Vector{Call} # keep track of queued calls. Calls can be queued after call arrivalTime + dispatchDelay
	
	resim::Resimulation
	
	mobilisationDelay::MobilisationDelay
	
	# decision logic
	addCallToQueue!::Function
	findAmbToDispatch!::Function
	redispatch::Redispatch
	moveUpData::MoveUpData
	
	# demand
	demand::Demand
	demandCoverage::DemandCoverage
	
	responseTravelPriorities::Dict{Priority,Priority} # responseTravelPriorities[p] gives the travel priority for responding to call of priority p
	targetResponseDurations::Vector{Float} # targetResponseDurations[Int(priority)] gives maximum desired response duration for call of given priority
	
	# for animation:
	currentCalls::Set{Call} # all calls between arrival and service finish at current time
	previousCalls::Set{Call} # calls in currentCalls for previous frame
	
	stats::SimStats
	
	reps::Vector{Simulation} # simulation replications, should have calls and stats defined
	numReps::Int # length(reps)
	isRunnable::Bool # whether or not simulation object can be run / simulated
	
	# files/folders:
	inputPath::String
	outputPath::String
	inputFiles::Dict{String,File} # given input file name (e.g. "ambulances"), returns file information
	outputFiles::Dict{String,File}
	
	eventsFile::EventsFile
	
	writeOutput::Bool # true if writing output files during simulation (e.g. events file), false otherwise (e.g. animating)
	initialised::Bool # true if simulation has been initialised and so can be run, false otherwise
	used::Bool # true if simulation has started running (and so restarting would require copying from backup)
	complete::Bool # true if simulation has ended (no events remaining)
	animating::Bool # true if being used for animation, false otherwise
	
	backup::Simulation # copy of simulation, for restarts (does not include a backup of all fields in order to save on memory, see backup! function for missing fields)
	
	configRootElt::XMLElement
	
	Simulation() = new(nullTime, nullTime, nullTime,
		Network(), Travel(), Map(), Grid(),
		[], [], [], [],
		0, 0, 0, 0,
		[], 0, [],
		Resimulation(),
		MobilisationDelay(),
		nullFunction, nullFunction, Redispatch(Val{:default}), MoveUpData(),
		Demand(), DemandCoverage(),
		Dict(), [],
		Set(), Set(),
		SimStats(),
		[], 0, true,
		"", "", Dict(), Dict(),
		EventsFile(),
		false, false, false, false, false)
end
