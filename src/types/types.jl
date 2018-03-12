type Location
	x::Float # latitude, or other
	y::Float # longitude, or other
	
	Location() = new(nullX, nullY)
	Location(x,y) = new(x,y)
end

type Node
	index::Int
	location::Location
	offRoadAccess::Bool # if node can be used to get on-road and off-road
	
	Node() = new(nullIndex, Location(), true)
end

type Arc
	index::Int
	fromNodeIndex::Int
	toNodeIndex::Int
	
	Arc() = new(nullIndex, nullIndex, nullIndex)
end

# graph data for network (actually a digraph)
type Graph
	# parameters:
	isReduced::Bool
	nodes::Vector{Node}
	arcs::Vector{Arc}
	
	light::LightGraphs.DiGraph
	fadjList::Vector{Vector{Int}} # shorthand for light.fadjlist; fadjList[i] gives nodes that are connected to node i by an outgoing arc from node i
	
	# for full graph:
	nodePairArcIndex::SparseMatrixCSC{Int,Int} # for node indices i,j, nodePairArcIndex[i,j] should return arc index for node pair
	# cannot always use nodePairArcIndex for reduced graph, there may be more than one arc from node i to j, use spNodePairArcIndex instead
	
	Graph(isReduced::Bool) = new(isReduced, [], [],
		LightGraphs.DiGraph(), [],
		spzeros(Int, 0, 0))
end

# travel data for network
type NetTravel
	# parameters:
	isReduced::Bool
	modeIndex::Int
	arcTimes::Vector{Float} # arcTimes[i] gives travel time along arc i
	
	# for use in reduced graph:
	spTimes::Array{FloatSpTime,2} # spTimes[i,j] = shortest path time between node i and j
	spFadjIndex::Array{IntFadj,2} # for shortest path from rNode i to j, spFadjIndex[i,j] gives the index (in fadjList[i], see Graph) of the successor rNode of i for this path
	spNodePairArcIndex::SparseMatrixCSC{Int,Int} # spNodePairArcIndex[i,j] = index of arc incident to nodes i,j, if it provides shortest travel time
	spFadjArcList::Vector{Vector{Int}} # spFadjArcList[i][j] gives index of rArc from rNode[i] to jth outoing rNode of rNode[i] (i.e. rGraph.fadjList[i][j])
	
	# for use in full graph:
	fNodeToRNodeTime::Vector{Dict{Int,Float}} # fNodeToRNodeTime[i][j] gives time from fNode[i] to rNode[j], as long as rNode[j] is in fNodeToRNodes[i] (see type Network)
	fNodeFromRNodeTime::Vector{Dict{Int,Float}} # fNodeFromRNodeTime[i][j] gives time to fNode[i] from rNode[j], as long as rNode[j] is in fNodeFromRNodes[i] (see type Network)
	rArcFNodesTimes::Vector{Vector{Float}} # rArcFNodesTimes[i][k] gives travel time from rArc[i].fromNodeIndex to rArcFNodes[i][k]
	
	# useful object to/from node data, for full graph (see commonFNodes in Network):
	commonFNodeToFNodeTime::Array{Float,2} # commonFNodeToFNodeTime[i,j] gives shortest path time from commonFNodes[i] to fNode j
	fNodeToCommonFNodeTime::Array{Float,2} # fNodeToCommonFNodeTime[i,j] gives shortest path time from fNode i to commonFNodes[j]
	commonFNodeToFNodeRNodes::Array{Vector{Int},2} # commonFNodeToFNodeRNodes[i,j] gives start and end rNode for shortest path from commonFNodes[i] to fNode j
	fNodeToCommonFNodeRNodes::Array{Vector{Int},2} # fNodeToCommonFNodeRNodes[i,j] gives start and end rNode for shortest path from fNode i to commonFNodes[j]
	fNodeNearestHospitalIndex::Vector{Int} # fNodeNearestHospitalIndex[i] gives index of nearest hospital from fNode[i]
	
	NetTravel(isReduced::Bool) = new(isReduced, nullIndex, [],
		Array{FloatSpTime,2}(0,0), Array{IntFadj,2}(0,0), spzeros(Int, 0, 0), [],
		[], [], [],
		Array{Float,2}(0,0), Array{Float,2}(0,0), Array{Vector{Int},2}(0,0), Array{Vector{Int},2}(0,0), [])
end

type Network
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
	rArcFNodeIndex::Vector{Dict{Int,Int}} # rArcFNodeIndex[i][j] gives index that fNode[j] appears in rArc[i]; should be same as find(rArcFNodes[i] .== j), so rArcFNodes[i][rArcFNodeIndex[i][j]] = j
	fNodeToRNodes::Vector{Vector{Int}} # fNodeToRNodes[i] gives indices of rNodes that fNode[i] can travel to, and is "near" ("near" meaning that the fNode is on an rArc incident to rNode, or fNode is rNode in rGraph; it is similar to adjacency)
	fNodeFromRNodes::Vector{Vector{Int}} # fNodeFromRNodes[i] gives indices of rNodes that can travel to fNode[i], and is "near"
	fNodeToRNodeNextFNode::Vector{Dict{Int,Int}} # fNodeToRNodeNextFNode[i][j] gives index of fNode after fNode[i] on path to rNode[j], where j is in fNodeToRNodes[i]. Needed to find path from an fNode to an rNode
	
	# fNodes that are common start/end points for travel (e.g, fNodes nearest to stations, hospitals):
	commonFNodes::Vector{Int} # list of common fNodes; = find(isFNodeCommon)
	isFNodeCommon::Vector{Bool} # isFNodeCommon[i] = true if fNode i is common, false otherwise
	fNodeCommonFNodeIndex::Vector{Int} # fNodeCommonFNodeIndex[i] gives index of fNode i in commonFNodes, if isFNodeCommon[i] = true
	
	Network() = new(Graph(false), Graph(true),
		[], [],
		[], [], [], [], [], [], [], [],
		[], [], [])
end

# travel information, for on-road and off-road
type TravelMode
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
type Travel
	numModes::Int # number of travel modes
	numSets::Int # number of sets of travel modes; each set may contain multiple travel modes e.g. for different combinations of travel priorities and ambulance classes. Different sets can overlap, by containing the same travel modes.
	
	modes::Vector{TravelMode}
	modeLookup::Array{Int,2} # modeLookup[i,j] = index of travel mode to use for travel set i, travel priority j. Change this variable according to modelling needs
	
	setsStartTimes::Vector{Float} # setsStartTimes[i] gives time at which travel set setsTimeOrder[i] should be started (setsStartTimes[i+1] or Inf gives end time)
	setsTimeOrder::Vector{Int} # setsTimeOrder[i] gives travel set to start using at time setsStartTimes[i]
	recentSetsStartTimesIndex::Int # index of most recently used value in setsStartTimes (and setsTimeOrder), should only ever increase in value
	
	Travel() = new(nullIndex, nullIndex,
		[], Array{Int,2}(0,0),
		[], [], nullIndex)
end

# for storing ambulance routes
# routes include a path (on fGraph), and may start/end off the graph
type Route
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
	startRNode::Int # index of first rNode in route
	startRNodeTime::Float # time at which startRNode is reached
	endRNode::Int # index of last rNode in route
	endRNodeTime::Float # time at which endRNode is reached
	endFNode::Int # index of last fNode in route
	endFNodeTime::Float # time at which endFNode is reached
	
	firstRArc::Int # index of first rArc in route; = nullIndex if startFNode == endFNode
	
	## fields that vary throughout the route are below here
	## most of the route code assumes that these fields only change to move "forward" through the route
	
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
	
	Route() = new(nullPriority, nullIndex,
		Location(), nullTime, Location(), nullTime,
		nullIndex, nullTime, nullIndex, nullTime, nullIndex, nullTime, nullIndex, nullTime,
		nullIndex,
		nullIndex, nullTime, nullTime,
		nullIndex, nullIndex, nullTime,
		nullIndex, nullIndex, nullTime)
end

type Event
	form::EventForm # "type" is taken
	time::Float
	ambIndex::Int
	callIndex::Int
	stationIndex::Int # for now, only use this for resimulation, otherwise use ambulances[ambIndex].stationIndex
	
	Event() = new(nullEvent, nullTime, nullIndex, nullIndex, nullIndex)
end

type Ambulance
	index::Int
	status::AmbStatus
	stationIndex::Int
	callIndex::Int
	route::Route
	event::Event # next/current event, useful for deleting next event from eventList
	class::AmbClass # class/type of ambulance
	# schedule::Schedule # to be added
	
	# for animation:
	currentLoc::Location
	movedLoc::Bool
	
	# for statistics:
	totalTravelTime::Float # this should only be updated after finishing each route / sim end
	totalBusyTime::Float # total time that ambulance has been busy
	# totalStationTime::Float # total time spent at station
	numCallsTreated::Int # total number of calls that ambulance provided treatment for
	numCallsTransferred::Int # total number of calls transferred to hospital
	numDiversions::Int # number of times that ambulance is diverted from one call to another
	atStationDispatches::Int # total number of dispatches while at station
	onRoadDispatches::Int # total number of dispatches while on road
	afterServiceDispatches::Int # total number of dispatches directly after providing service at callout
	
	Ambulance() = new(nullIndex, ambNullStatus, nullIndex, nullIndex, Route(), Event(), nullAmbClass,
		Location(), false,
		0.0, 0.0, 0, 0, 0, 0, 0, 0)
end

type Call
	index::Int
	status::CallStatus
	ambIndex::Int
	priority::Priority
	transfer::Bool # true if requires transfer to hospital
	hospitalIndex::Int # hospital (if any) that ambulance is transferred to. If hospitalIndex == nullIndex, will transfer to nearest hospital
	location::Location # where call occurs
	
	arrivalTime::Float # time at which call arrives
	dispatchDelay::Float # delay between call arrival and considering dispatching an ambulance
	onSceneDuration::Float # time spent at call location
	transferDuration::Float # for hospital transfer
	
	# node nearest to call location:
	nearestNodeIndex::Int
	nearestNodeDist::Float
	
	# for animation:
	currentLoc::Location
	movedLoc::Bool
	
	# for statistics:
	dispatchTime::Float # time at which final ambulance was dispatched (= arrivalTime + dispatchDelay, unless call queued/bumped)
	ambArrivalTime::Float # time at which ambulance arrives on-site
	responseTime::Float # time (duration) between call arrival and ambulance arrival at call location
	hospitalArrivalTime::Float # time at which ambulance arrives at hospital
	numBumps::Int # total number of times that call gets bumped
	wasQueued::Bool # whether call was queued or not
	
	Call() = new(nullIndex, callNullStatus, nullIndex, nullPriority, true, nullIndex, Location(),
		nullTime, nullTime, nullTime, nullTime,
		nullIndex, nullDist,
		Location(), false,
		nullTime, nullTime, nullTime, nullTime, 0, false)
end

type Hospital
	index::Int
	location::Location
	nearestNodeIndex::Int
	nearestNodeDist::Float
	
	# for statistics:
	numTransfers::Int # total number of patient transfers from ambulances to hospital
	
	Hospital() = new(nullIndex, Location(), nullIndex, nullDist,
		0)
end

type Station
	index::Int
	location::Location
	capacity::Int # maximum number of ambulances that station can hold
	nearestNodeIndex::Int
	nearestNodeDist::Float
	
	# for statistics:
	# totalAmbIdleTime::Float # total time that ambulances are idle at station
	
	Station() = new(nullIndex, Location(), 0, nullIndex, nullDist)
end

type Map
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
type GridSearchRect
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
type GridRect
	nodeIndices::Vector{Int}
	
	GridRect() = new([])
end

# grid breaks map region into rectangles,
# each rectangle stores graph node indices (from full graph)
# this is for quickly finding the nearest node to a location
type Grid
	nx::Int # number of divisions in x direction
	ny::Int # number of divisions in y direction
	xRange::Float # size of divisions in x direction, = map.xRange / nx
	yRange::Float # size of divisions in y direction, = map.yRange / ny
	xRangeDist::Float # distance of divisions in x direction, = xRange * map.xScale
	yRangeDist::Float # distance of divisions in y direction, = yRange * map.yScale
	rects::Array{GridRect,2} # rectangles
	searchRect::GridSearchRect
	
	Grid() = new(nullIndex, nullIndex, nullDist, nullDist, nullDist, nullDist, Array{GridRect,2}(0, 0), GridSearchRect())
	function Grid(map::Map, nx, ny)
		grid = new(nx, ny, nullDist, nullDist, nullDist, nullDist, Array{GridRect,2}(nx, ny))
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
type Raster
	# parameters:
	x::Vector{Float}
	y::Vector{Float}
	z::Array{Float,2} # z[i,j] corresponds with x[i], y[j]
	
	nx::Int # length(x)
	ny::Int # length(y)
	dx::Float # x step size
	dy::Float # y step size
	
	function Raster(x, y, z)
		nx = length(x)
		ny = length(y)
		assert(nx > 1)
		assert(ny > 1)
		assert((nx, ny) == size(z))
		dx = (maximum(x) - minimum(x))/(nx-1)
		dy = (maximum(y) - minimum(y))/(ny-1)
		assert(dx > 0)
		assert(dy > 0)
		return new(x, y, z, nx, ny, dx, dy)
	end
end

# move up data types
abstract type MoveUpDataType end
type EmptyMoveUpData <: MoveUpDataType end

# compliance table data
type CompTableData <: MoveUpDataType
	# parameters:
	compTable::Array{Int,2} # compTable[i,j] = number of ambulances to place at station j, with i total idle ambs
	
	# arrays for recycling:
	ambMovable::Vector{Bool} # ambMovable[i] = true if ambulance i is available for move up, false otherwise
	
	CompTableData() = new(Array{Int,2}(0,0),
		[])
end

# dmexclp - dynamic maximum expected coverage location problem
type DmexclpData <: MoveUpDataType
	# parameters:
	coverTime::Float # ambulance 'covers' a location if it can reach the location in this time
	busyFraction::Float # fraction for which each ambulance is busy, approximate
	
	marginalBenefit::Vector{Float} # marginalBenefit[i] = benefit of adding an ith ambulance to cover single demand
	
	# arrays for recycling:
	nodeDemand::Vector{Int} # nodeDemand[i] = demand at node i
	stationCoverNode::Array{Bool,2} # stationCoverNode[i,j] = true if station i 'covers' node j
	stationCoverNodes::Vector{Vector{Int}} # stationCoverNodes[i] = list of node indices covered by station i
	stationNumIdleAmbs::Vector{Int} # number of idle ambulances assigned to each station
	nodeCoverCount::Vector{Int} # nodeCoverCount[i] = number of idle ambulances covering node i
	
	DmexclpData() = new(nullTime, 0.0,
		[],
		[], Array{Bool,2}(0,0), [], [], [])
end

type PriorityListData <: MoveUpDataType
	# parameters:
	priorityList::Vector{Int} # priorityList[i] gives station index that the ith free ambulance should be moved to
	
	# arrays for recycling:
	stationNumIdleAmbs::Vector{Int} # number of idle ambulances assigned to each station
	
	PriorityListData() = new([],
		[])
end

type ZhangIpData <: MoveUpDataType
	# parameters:
	busyFraction::Float # ambulance busy fraction - should remove this, make marginalBenefit a parameter
	travelTimeCost::Float # travel time cost multiplier
	maxIdleAmbTravelTime::Float # max travel time for idle ambulances. 0.021 (days) is about 30 minutes
	maxNumNearestStations::Int # number of nearest stations to consider for each ambulance (may include ambulance's station)
	
	marginalBenefit::Vector{Float} # marginalBenefit[i] = benefit of adding an ith ambulance to a station
	
	ZhangIpData() = new(0.0, 0.0, Inf, 0,
		[])
end

type Temp1Data <: MoveUpDataType
	# parameters:
	benefit::Array{Vector{Float},2} # benefit[i1,i2][k] = benefit from having k ambulances collectively at stations i1 and i2
	busyFraction::Float # ambulance busy fraction - should remove this, make marginalBenefit a parameter
	stationPairs::Vector{Vector{Int}} # stationPairs[i] gives ith station pair, which is a pair of station indices
	travelTimeCost::Float # travel time cost multiplier
	maxIdleAmbTravelTime::Float # max travel time for idle ambulances
	maxNumNearestStations::Int # number of nearest stations to consider for each ambulance (may include ambulance's station)
	
	marginalBenefit::Array{Vector{Float},2} # marginal benefit values from ambulance to station allocations, calculated from 'benefit'
	
	Temp1Data() = new(Array{Vector{Float},2}(0,0), 0.0, Vector{Vector{Int}}(0), 0.0, Inf, 0,
		Array{Vector{Float},2}(0,0))
end

type Temp2Data <: MoveUpDataType
	# parameters:
	benefit::Array{Array{Float,2},2} # benefit[i1,i2][k1,k2] = benefit from having k1 and k2 ambulances at stations i1 and i2, respectively
	busyFraction::Float # ambulance busy fraction - should remove this, make marginalBenefit a parameter
	stationPairs::Vector{Vector{Int}} # stationPairs[i] gives ith station pair, which is a pair of station indices
	travelTimeCost::Float # travel time cost multiplier
	maxIdleAmbTravelTime::Float # max travel time for idle ambulances
	maxNumNearestStations::Int # number of nearest stations to consider for each ambulance (may include ambulance's station)
	
	marginalBenefit::Array{Array{Float,2},2} # marginal benefit values from ambulance to station allocations, calculated from 'benefit'
	
	Temp2Data() = new(Array{Array{Float,2},2}(0,0), 0.0, Vector{Vector{Int}}(0), 0.0, Inf, 0,
		Array{Array{Float,2},2}(0,0))
end

type MoveUpData
	# parameters:
	useMoveUp::Bool
	moveUpModule::MoveUpModule # indicates move up module to be used
	
	# move up data types:
	compTableData::CompTableData
	dmexclpData::DmexclpData
	priorityListData::PriorityListData
	zhangIpData::ZhangIpData
	temp1Data::Temp1Data
	temp2Data::Temp2Data
	
	MoveUpData() = new(false, nullMoveUpModule,
		CompTableData(), DmexclpData(), PriorityListData(), ZhangIpData(), Temp1Data(), Temp2Data())
end

type File
	name::String
	path::String # includes name
	iostream::IOStream
	checksum::UInt
	
	# File() = new()
	File() = new("", "", IOStream(""), 0)
end

type Resimulation
	# parameters:
	use::Bool # true if resimulating (will follow event trace), false otherwise
	timeTolerance::Float
	
	events::Vector{Event}
	prevEventIndex::Int # index of previous event in events field
	
	Resimulation() = new(false, 0.0,
		[], nullIndex)
end

type Simulation
	startTime::Float
	time::Float
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
	
	eventList::Vector{Event}
	queuedCallList::Vector{Call} # keep track of queued calls. Calls can be queued after call arrivalTime + dispatchDelay
	
	resim::Resimulation
	
	# for move up:
	moveUpData::MoveUpData
	
	# for statistics:
	targetResponseTimes::Vector{Float} # targetResponseTimes[Int(priority)] gives maximum desired response time for call of given priority
	
	# for animation:
	currentCallList::Vector{Call} # all calls between arrival and service finish at current time
	previousCallList::Vector{Call} # calls in currentCallList for previous frame
	
	# files/folders:
	inputPath::String
	outputPath::String
	inputFiles::Dict{String,File} # given input file name (e.g. "ambulances"), returns file information
	outputFiles::Dict{String,File}
	eventsFileIO::IOStream # write simulation trace of events to this file
	
	writeOutput::Bool # true if outputFiles should be used to write output (false if animating)
	
	used::Bool # true if simulation has started running (and so restarting would require copying from backup)
	complete::Bool # true if simulation has ended (no events remaining)
	backup::Simulation # copy of simulation, for restarts (does not include net, travel, grid, or resim, as copying these would waste memory)
	
	configRootElt::XMLElement
	
	Simulation() = new(nullTime, nullTime, nullTime,
		Network(), Travel(), Map(), Grid(),
		[], [], [], [],
		[], [],
		Resimulation(),
		MoveUpData(),
		[],
		[], [],
		"", "", Dict(), Dict(), IOStream(""),
		false,
		false, false)
end
