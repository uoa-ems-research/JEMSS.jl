# initialise travel data for fGraph
# arcTravelTimes[i,j] gives travel time along arc j (in fGraph) for travel mode i
function initFNetTravels!(net::Network, arcTravelTimes::Array{Float,2})
	numModes = size(arcTravelTimes, 1) # shorthand
	assert(length(net.fGraph.arcs) == size(arcTravelTimes, 2)) # check number of arcs
	
	# create fNetTravels
	net.fNetTravels = [NetTravel(false) for i = 1:numModes]
	for i = 1:numModes
		fNetTravel = net.fNetTravels[i]
		fNetTravel.modeIndex = i
		fNetTravel.arcTimes = arcTravelTimes[i,:]
	end
end

# return true if fNode i is in rGraph, false otherwise
function isFNodeInRGraph(net::Network, i::Int)
	return (net.fNodeRNode[i] != nullIndex)
end

# create the reduced graph from the full graph, for the given network
function createRGraphFromFGraph!(net::Network)
	# assumptions:
	# - no nodes are connected to themselves
	# - there is at least one intersection node (node adjacent to more than 2 other nodes)
	
	net.rGraph = Graph(true)
	
	# shorthand names:
	fGraph = net.fGraph
	fNodes = fGraph.nodes
	numFNodes = length(fNodes)
	rGraph = net.rGraph
	
	##############################################################################
	
	# some useful fGraph connectivity data
	
	# shorthand:
	fNodesOut = fGraph.light.fadjlist # fNodesOut[i] gives indices of fNodes incident to arcs outgoing from fNode i
	fNodesIn = fGraph.light.badjlist # fNodesIn[i] gives indices of fNodes incident to arcs incoming to fNode i
	
	fNodeNumAdjacent = [length(union(fNodesIn[i], fNodesOut[i])) for i = 1:numFNodes] # number of fNodes adjacent to each fNode
	
	# check that no node is connected to itself
	for i = 1:numFNodes
		assert(!in(i, fNodesIn[i]))
	end
	
	##############################################################################
	
	# find intersection/decision nodes in fGraph, add to rGraph
	
	# only keep nodes with >= 3 adjacent nodes, or 2 adjacent nodes if the incoming and outgoing arcs cannot be replaced
	# Also, for easier programming, put leaf nodes in fGraph into rGraph, hopefully there are not too many...
	keepNode = Vector{Bool}(numFNodes) # find which fNodes belong in rGraph
	keepNode[:] = false
	for i = 1:numFNodes
		if fNodeNumAdjacent[i] >= 3
			keepNode[i] = true
		elseif fNodeNumAdjacent[i] == 2 && (length(fNodesIn[i]) + length(fNodesOut[i]) == 3)
			keepNode[i] = true
		elseif fNodeNumAdjacent[i] == 1
			keepNode[i] = true
		end
	end
	rGraph.nodes = deepcopy(fNodes[keepNode]) # need to renumber rNode indices later
	net.rNodeFNode = find(keepNode)
	net.fNodeRNode = [nullIndex for i = 1:numFNodes]
	net.fNodeRNode[keepNode] = [1:sum(keepNode);]
	
	# special consideration: need to add one decision node to any loop with currently only one rNode
	# Otherwise, fNodes in the loop still have a decision as to how to reach the rNode
	for ri = 1:length(rGraph.nodes)
		fi = net.rNodeFNode[ri]
		for fj = fNodesOut[fi] # fNodes outgoing from fNode[fi]
			firstFNode = fj
			prevFNode = fi
			fNode = fj
			# traverse fGraph from rNode[ri], until another rNode or dead end (leaf node?) reached
			while !isFNodeInRGraph(net, fNode) # !(isFNodeInRGraph(net, fNode) || fNodeNumAdjacent[fNode] == 1)
				# traverse along an arc to next node
				nextFNode = first(setdiff(fNodesOut[fNode], prevFNode)) # fNodesOut[fNode] = indices of fNodes outgoing from fNode[fNode]
				prevFNode = fNode
				fNode = nextFNode
			end
			# if we have looped back to same rNode, add firstFNode as decision nodes
			if ri == net.fNodeRNode[fNode]
				keepNode[firstFNode] = true
			end
		end
	end
	
	rGraph.nodes = deepcopy(fNodes[keepNode])
	# renumber rNode indices
	for i = 1:length(rGraph.nodes)
		rGraph.nodes[i].index = i
	end
	net.rNodeFNode = find(keepNode)
	net.fNodeRNode = [nullIndex for i = 1:numFNodes]
	net.fNodeRNode[keepNode] = [1:sum(keepNode);]
	
	# shorthand names:
	rNodes = rGraph.nodes
	fNodeRNode = net.fNodeRNode
	rNodeFNode = net.rNodeFNode
	numRNodes = length(rNodes)
	
	assert(numRNodes > 0) # may need a different createRGraphFromFGraph!() function if numRNodes == 0
	
	# check index conversion between fNodes and rNodes
	assert(all([i == fNodeRNode[rNodeFNode[rNodes[i].index]] for i = 1:numRNodes]))
	# assert(all([i == fNodeRNode[rNodeFNode[i]] for i = 1:numRNodes]))
	
	##############################################################################
	
	# calculate rArcs, rArcFNodes, fNodeRArcs
	# find outgoing arcs from each rNode to other rNodes
	net.fNodeFromRNodes = [Vector{Int}(0) for i = 1:numFNodes]
	net.rArcFNodes = []
	net.fNodeRArcs = [Vector{Int}(0) for i = 1:numFNodes]
	rArcFNodes = net.rArcFNodes # shorthand
	fNodeRArcs = net.fNodeRArcs # shorthand
	rArcs = rGraph.arcs # shorthand
	numRArcs = 0
	for ri = 1:numRNodes
		fi = rNodeFNode[ri]
		for fj = fNodesOut[fi] # fNodes outgoing from fNode[fi]
			numRArcs += 1 # will create another rArc
			# traverse fGraph from rNode[ri], until another rNode or dead end (leaf node?) reached
			fNodesOnRArc = [fi] # indices of fNodes on rArc
			prevFNode = fi
			fNode = fj
			while !isFNodeInRGraph(net, fNode)
				push!(fNodeRArcs[fNode], numRArcs)
				push!(fNodesOnRArc, fNode)
				# traverse along an arc to next node
				nextFNode = first(setdiff(fNodesOut[fNode], prevFNode)) # fNodesOut[fNode] = indices of fNodes outgoing from fNode
				prevFNode = fNode
				fNode = nextFNode
			end
			push!(fNodesOnRArc, fNode)
			rj = fNodeRNode[fNode] # rNode at end of rArc, fNodeRNode[fNode] = nullIndex if fNode[fNode] is not in rGraph
			assert(rj != nullIndex) # should not happen; earlier we put all leaf nodes in fGraph in rGraph
			assert(ri != rj) # ri == rj would mean that arc is a loop, should have added enough rNodes to remove this problem earlier
			
			# add rArc
			arc = Arc()
			arc.index = numRArcs
			arc.fromNodeIndex = ri
			arc.toNodeIndex = rj
			push!(rArcs, arc)
			push!(rArcFNodes, fNodesOnRArc)
		end
	end
	
	assert(numRArcs > 0)
	
	# calculate rArcFNodeIndex from rArcFNodes
	net.rArcFNodeIndex = [Dict{Int,Int}() for i = 1:numRArcs]
	for ri = 1:numRArcs
		for (i, fi) in enumerate(rArcFNodes[ri])
			net.rArcFNodeIndex[ri][fi] = i
		end
	end
	
	# set fNodeRArcs for each rNode (for corresponding fNode)
	for i = 1:numRArcs
		fromRNode = rArcs[i].fromNodeIndex
		toRNode = rArcs[i].toNodeIndex
		push!(fNodeRArcs[rNodeFNode[fromRNode]], i)
		push!(fNodeRArcs[rNodeFNode[toRNode]], i)
	end
	
	# check that fNodeRArcs is not empty for rNodes, and other nodes have 1-2 entries
	for fi = 1:numFNodes
		if isFNodeInRGraph(net, fi)
			assert(!isempty(fNodeRArcs[fi]))
		else
			assert(length(fNodeRArcs[fi]) >= 1 && length(fNodeRArcs[fi]) <= 2)
		end
	end
	
	# check that start and end fNodes on each rArc are also in rGraph, and arcs are numbered 1 to n
	for i = 1:numRArcs
		assert(isFNodeInRGraph(net, rArcFNodes[i][1]))
		assert(isFNodeInRGraph(net, rArcFNodes[i][end]))
		assert(rArcs[i].index == i)
	end
	
	##############################################################################
	
	# calculate fNodeToRNodes and fNodeFromRNodes
	
	net.fNodeToRNodes = [Vector{Int}(0) for i = 1:numFNodes]
	net.fNodeFromRNodes = [Vector{Int}(0) for i = 1:numFNodes]
	net.fNodeToRNodeNextFNode = [Dict{Int,Int}() for i = 1:numFNodes]
	
	# shorthand:
	fNodeToRNodes = net.fNodeToRNodes
	fNodeFromRNodes = net.fNodeFromRNodes
	fNodeToRNodeNextFNode = net.fNodeToRNodeNextFNode
	
	# traverse each rArc
	for rArc = rArcs
		fNodesOnRArc = rArcFNodes[rArc.index]
		assert(isFNodeInRGraph(net, fNodesOnRArc[1]))
		assert(isFNodeInRGraph(net, fNodesOnRArc[end]))
		startRNode = fNodeRNode[fNodesOnRArc[1]]
		endRNode = fNodeRNode[fNodesOnRArc[end]]
		# for fNodes on arc (excluding rNodes):
		for i = 2:length(fNodesOnRArc)-1
			fNode = fNodesOnRArc[i]
			push!(fNodeToRNodes[fNode], endRNode)
			push!(fNodeFromRNodes[fNode], startRNode)
			fNodeToRNodeNextFNode[fNode][endRNode] = fNodesOnRArc[i+1]
		end
	end
	# for rNodes, set corresponding fNode to have fNodeToRNodes and fNodeFromRNodes return the rNode
	for ri = 1:numRNodes
		fi = rNodeFNode[ri]
		fNodeFromRNodes[fi] = [ri]
		fNodeToRNodes[fi] = [ri]
	end
	
	# check that fNodeFromRNodes and fNodeToRNodes have been filled
	for i = 1:numFNodes
		assert(!isempty(fNodeFromRNodes[i]))
		assert(!isempty(fNodeToRNodes[i]))
	end
	
	##############################################################################
	
	# calculate fields for rGraph: light
	initGraph!(rGraph)
end

# calculate arcTimes, spTimes, and spSuccs, for rNetTravels
# also, calculate fNodeToRNodeTime, and fNodeFromRNodeTime, for fNetTravels...
# can only call after rGraph is created
function createRNetTravelsFromFNetTravels!(net::Network)
	
	# shorthand names:
	fGraph = net.fGraph
	rGraph = net.rGraph
	fNetTravels = net.fNetTravels
	numTravelModes = length(fNetTravels)
	rArcFNodes = net.rArcFNodes
	numFNodes = length(fGraph.nodes)
	numRNodes = length(rGraph.nodes)
	numRArcs = length(rGraph.arcs)
	fNodeRNode = net.fNodeRNode
	rNodeFNode = net.rNodeFNode
	fNodeToRNodes = net.fNodeToRNodes
	
	rNetTravels = net.rNetTravels = [NetTravel(true) for i = 1:numTravelModes]
	
	assert(length(fGraph.arcs) == length(fNetTravels[1].arcTimes))
	
	for travelModeIndex = 1:numTravelModes
		# shorthand:
		fNetTravel = fNetTravels[travelModeIndex]
		rNetTravel = rNetTravels[travelModeIndex]
		
		rNetTravel.arcTimes = Vector{Float}(numRArcs)
		fNetTravel.fNodeToRNodeTime = [Dict{Int, Float}() for i = 1:numFNodes]
		fNetTravel.fNodeFromRNodeTime = [Dict{Int, Float}() for i = 1:numFNodes]
		fNetTravel.rArcFNodesTimes = [Vector{Float}() for i = 1:numRArcs]
		fNodeToRNodeTime = fNetTravel.fNodeToRNodeTime # shorthand
		fNodeFromRNodeTime = fNetTravel.fNodeFromRNodeTime # shorthand
		rArcFNodesTimes = fNetTravel.rArcFNodesTimes
		
		rNetTravel.modeIndex = fNetTravel.modeIndex
		
		for rArcIndex = 1:numRArcs
			fNodesOnRArc = rArcFNodes[rArcIndex]
			startRNode = fNodeRNode[fNodesOnRArc[1]]
			endRNode = fNodeRNode[fNodesOnRArc[end]]
			
			# traverse arc, calculate rArcTime and fNodeFromRNodeTime
			rArcTime = 0.0
			for i = 1:length(fNodesOnRArc)-1
				fNode = fNodesOnRArc[i]
				nextFNode = fNodesOnRArc[i+1]
				fNodeFromRNodeTime[fNode][startRNode] = rArcTime
				push!(rArcFNodesTimes[rArcIndex], rArcTime)
				fArcIndex = fGraph.nodePairArcIndex[fNode, nextFNode]
				rArcTime += fNetTravel.arcTimes[fArcIndex]
			end
			push!(rArcFNodesTimes[rArcIndex], rArcTime)
			rNetTravel.arcTimes[rArcIndex] = rArcTime
			
			# calculate fNodeToRNodeTime
			for i = 2:length(fNodesOnRArc)-1
				fNode = fNodesOnRArc[i]
				fNodeToRNodeTime[fNode][endRNode] = rArcTime - fNodeFromRNodeTime[fNode][startRNode]
			end
			fNodeToRNodeTime[fNodesOnRArc[end]][endRNode] = 0.0 # time to self is 0
			
			# check length of rArcFNodesTimes is same as length of rArcFNodes for current arc
			assert(length(fNodesOnRArc) == length(rArcFNodesTimes[rArcIndex]))
		end
		
		# # for rNodes, set corresponding fNode to have fNodeFromRNodeTime and fNodeToRNodeTime return 0
		# for ri = 1:numRNodes
			# fi = rNodeFNode[ri]
			# # fNodeToRNodeTime[fi][ri] = 0.0 # should already be true?
			# # fNodeFromRNodeTime[fi][ri] = 0.0 # should already be true?
		# end
		
		# for rNodes, check that corresponding fNode has fNodeFromRNodeTime and fNodeToRNodeTime return 0
		for ri = 1:numRNodes
			fi = rNodeFNode[ri]
			assert(fNodeToRNodeTime[fi][ri] == 0.0)
			assert(fNodeFromRNodeTime[fi][ri] == 0.0)
		end
		
		# for fNodes along each rArc, fNodeFromRNodeTime and fNodeToRNodeTime should add to rArc time
		for rArcIndex = 1:numRArcs
			rArcTime = rNetTravel.arcTimes[rArcIndex]
			fNodesOnRArc = rArcFNodes[rArcIndex]
			startRNode = fNodeRNode[fNodesOnRArc[1]]
			endRNode = fNodeRNode[fNodesOnRArc[end]]
			for i = 2:length(fNodesOnRArc)-1
				fNode = fNodesOnRArc[i]
				assert(abs(fNodeFromRNodeTime[fNode][startRNode] + fNodeToRNodeTime[fNode][endRNode] - rArcTime) < eps())
			end
		end
	end
	
	calcRNetTravelsShortestPaths!(net)
end

# For all rNetTravels, calculate and store shortest path data
# Uses repeated Dijkstra's shortest path algorithm,
# for sparse graphs, this is faster than Floyd-Warshall algorithm
function calcRNetTravelsShortestPaths!(net::Network)
	# shorthand
	rGraph = net.rGraph
	numRNodes = length(rGraph.nodes)
	rNetTravels = net.rNetTravels
	numTravelModes = length(rNetTravels)
	
	# initialise sp data
	for rNetTravel in rNetTravels
		rNetTravel.spTimes = Array{FloatSpTime,2}(numRNodes, numRNodes)
		rNetTravel.spFadjIndexDelta = [spzeros(IntFadj, numRNodes) for i = 1:numRNodes]
	end
	spFadjIndexBase = net.spFadjIndexBase = Array{IntFadj,2}(numRNodes, numRNodes)
	
	# create reverse graph
	revRGraph = LightGraphs.reverse(rGraph.light)
	
	# for each travel mode, create arc time matrix for reverse graph, for use by LightGraphs
	revArcTimeArrays = [spzeros(FloatSpTime, numRNodes, numRNodes) for i = 1:numTravelModes]
	for rNetTravel in rNetTravels
		revArcTimeArray = revArcTimeArrays[rNetTravel.modeIndex]
		for rArc in rGraph.arcs
			revArcTimeArray[rArc.toNodeIndex, rArc.fromNodeIndex] = rNetTravel.arcTimes[rArc.index]
		end
	end
	
	# calculate shortest paths for each possible destination node, within each travel mode
	# only consider one destination at a time, to reduce memory usage
	rNetTravelsSpFadjIndex = [Vector{IntFadj}(numRNodes) for i = 1:numTravelModes] # rNetTravelsSpFadjIndex[i][j] will store successor of node j (in fadjList[j]) to the current end node, for rNetTravels[i]
	fadjList = rGraph.fadjList # shorthand
	maxFadjLength = maximum([length(fadjList[i]) for i = 1:numRNodes]) # maximum number of arcs outgoing from any node in rGraph
	fadjIndexFreq = Vector{Int}(maxFadjLength) # fadjIndexFreq[i] will count frequency with which rNetTravelsSpFadjIndex[j][k] is i, for all travel modes
	for j = 1:numRNodes
		for rNetTravel in rNetTravels
			spData = LightGraphs.dijkstra_shortest_paths(revRGraph, j, revArcTimeArrays[rNetTravel.modeIndex])
			rNetTravel.spTimes[:,j] = spData.dists
			# spData.parents[i] is the node index to travel to after node i on path from i to j,
			# need to convert this to the index of the node in fadjList[i]
			for i = 1:numRNodes
				if i != j
					rNetTravelsSpFadjIndex[rNetTravel.modeIndex][i] = findfirst(fadjList[i], spData.parents[i])
				end
			end
		end
		
		# for each start node, find most common successor node over all rNetTravels,
		# (looking at rNetTravelsSpFadjIndex), and place in net.spFadjIndexBase;
		# for each travel mode with a different successor node, store in rNetTravel.spFadjIndexDelta
		for i = 1:numRNodes
			if i == j
				spFadjIndexBase[i,j] = nullIndex
			else
				# find most common fadjIndex among all travel modes
				fill!(fadjIndexFreq, 0)
				for rNetTravel in rNetTravels
					fadjIndexFreq[rNetTravelsSpFadjIndex[rNetTravel.modeIndex][i]] += 1
				end
				
				# set spFadjIndexBase[i,j] to store most common successor for path from node i to j
				fadjIndexBase = findmax(fadjIndexFreq)[2]
				spFadjIndexBase[i,j] = fadjIndexBase
				
				# for travel modes with a different successor node, set spFadjIndexDelta[i,j]
				for rNetTravel in rNetTravels
					fadjIndex = rNetTravelsSpFadjIndex[rNetTravel.modeIndex][i]
					if fadjIndex != fadjIndexBase
						rNetTravel.spFadjIndexDelta[i][j] = fadjIndex
					end
				end
			end
		end
	end
	
	# change sp successor data where needed
	# where a path i -> l has spNextNode(i,l) = k, but spNextNode(i,k) = j,
	# need to change all cases where spNextNode(i,-) = k to instead have spNextNode(i,-) = j
	# this case can happen when an arc (i,k) has the same travel time as a path (i,j,k)
	for rNetTravel in rNetTravels
		travelModeIndex = rNetTravel.modeIndex # shorthand
		for rArc in rGraph.arcs
			i = rArc.fromNodeIndex
			k = l = rArc.toNodeIndex
			j = shortestPathNextRNode(net, travelModeIndex, i, k)
			while j != k
				k = j
				j = shortestPathNextRNode(net, travelModeIndex, i, k)
			end
			if k != l
				# for any paths from i that have successor node l, change successor to be k
				for j = [1:i-1;i+1:numRNodes]
					if shortestPathNextRNode(net, travelModeIndex, i, j) == l
						setRNetTravelSpFadjIndexDelta!(net, rNetTravel, i, j, IntFadj(findfirst(fadjList[i], k)))
					end
				end
			end
		end
	end
	
	# set spNodePairArcIndex values
	for rNetTravel in rNetTravels
		travelModeIndex = rNetTravel.modeIndex # shorthand
		rNetTravel.spNodePairArcIndex = spzeros(Int, numRNodes, numRNodes)
		rNetTravel.spFadjArcList = deepcopy(rGraph.fadjList)
		for rArc in rGraph.arcs
			i = rArc.fromNodeIndex
			j = rArc.toNodeIndex
			if shortestPathNextRNode(net, travelModeIndex, i, j) == j &&
				rNetTravel.spTimes[i,j] == FloatSpTime(rNetTravel.arcTimes[rArc.index])
				
				rNetTravel.spNodePairArcIndex[i,j] = rArc.index
				rNetTravel.spFadjArcList[i][findfirst(rGraph.fadjList[i], j)] = rArc.index
			end
		end
	end
	
	# check stored shortest path times are same as those from traversing stored shortest paths
	# this can be slow, and should not be necessary unless shortest path code above has changed
	if checkMode && false
		for rNetTravel in rNetTravels
			checkShortestPathTimes(net, rNetTravel)
		end
	end
	
	# check assumptions
	for rNetTravel in rNetTravels
		for i = 1:numRNodes, j = [1:i-1;i+1:numRNodes]
			assert(0 < rNetTravel.spTimes[i,j] < Inf)
		end
	end
end

# For each rNetTravel, rebalance net.spFadjIndexBase, and rNetTravel.spFadjIndexDelta
# so that memory used to store shortest path data is minimised
# Only need to call this if storage of sp data was not minimised during calculation
function rebalanceRNetTravelsSpFadjIndices!(net::Network)
	# shorthand:
	rGraph = net.rGraph
	numRNodes = length(rGraph.nodes)
	rNetTravels = net.rNetTravels
	
	maxFadjLength = maximum([length(rGraph.fadjList[i]) for i = 1:numRNodes]) # maximum number of arcs outgoing from any node in rGraph
	fadjIndices = Vector{IntFadj}(maxFadjLength) # fadjIndices[i] will store the fadjIndex for rNetTravel i
	fadjIndexFreq = Vector{IntFadj}(maxFadjLength) # fadjIndexFreq[i] will count the frequency with which rNetTravel.spFadjIndexDelta[j][k] || net.spFadjIndexBase[j,k] is i, for all travel modes
	for i = 1:numRNodes, j = [1:i-1;i+1:numRNodes]
		# find the most common fadjIndex among all rNetTravels
		fill!(fadjIndexFreq, 0)
		for (k, rNetTravel) in enumerate(rNetTravels)
			fadjIndex = rNetTravel.spFadjIndexDelta[i][j]
			if fadjIndex == 0
				fadjIndex = net.spFadjIndexBase[i,j]
			end
			fadjIndices[k] = fadjIndex
			fadjIndexFreq[fadjIndex] += 1
		end
		
		# update fadjIndex values if most common index changed
		fadjIndexBase = findmax(fadjIndexFreq)[2]
		if net.spFadjIndexBase[i,j] != fadjIndexBase
			net.spFadjIndexBase[i,j] = fadjIndexBase
			
			# udpate spFadjIndexDelta[i][j] for each rNetTravel
			for (k, rNetTravel) in enumerate(rNetTravels)
				setRNetTravelSpFadjIndexDelta!(net, rNetTravel, i, j, fadjIndices[k])
				## or:
				# fadjIndex = fadjIndices[k]
				# if fadjIndex == fadjIndexBase
					# rNetTravel.spFadjIndexDelta[i][j] = 0
				# else
					# rNetTravel.spFadjIndexDelta[i][j] = fadjIndex
				# end
			end
		end
	end
	
	# # remove zero entries from sparse arrays
	# for rNetTravel in rNetTravels
		# for i = 1:numRNodes
			# dropzeros!(rNetTravel.spFadjIndexDelta[i])
		# end
	# end
end

function setRNetTravelSpFadjIndexDelta!(net::Network, rNetTravel::NetTravel, i::Int, j::Int, spFadjIndex::IntFadj)
	# assert(rNetTravel.isReduced)
	if net.spFadjIndexBase[i,j] == spFadjIndex && rNetTravel.spFadjIndexDelta[i][j] != 0
		rNetTravel.spFadjIndexDelta[i][j] = 0
		dropzeros!(rNetTravel.spFadjIndexDelta[i]) # not sure if this is slow
	else
		rNetTravel.spFadjIndexDelta[i][j] = spFadjIndex
	end
end

# For a given travel mode, check that the spSuccs
# gives paths with the same travel time as stored in spTimes
# For a path from node i to j, spSuccs[i,j] gives successor node index of i
function checkShortestPathTimes(rNetTravel::NetTravel, spSuccs::Array{IntRNode,2})
	# shorthand:
	spTimes = rNetTravel.spTimes
	numRNodes = size(spSuccs,1)
	
	assert(size(spSuccs) == (numRNodes, numRNodes))
	assert(size(rNetTravel.spTimes) == (numRNodes, numRNodes))
	
	# check each pair of start and end nodes
	for startRNode = 1:numRNodes, endRNode = 1:numRNodes
		if startRNode == endRNode
			assert(spTimes[startRNode, endRNode] == 0.0)
			assert(spSuccs[startRNode, endRNode] == nullIndex)
		else
			# follow path from start to end node, calculate and compare travel time
			i = startRNode
			t = 0.0
			while i != endRNode
				j = spSuccs[i,endRNode]
				rArcIndex = rNetTravel.spNodePairArcIndex[i,j]
				t += rNetTravel.arcTimes[rArcIndex]
				i = j # go to next node
			end
			assert(abs(spTimes[startRNode, endRNode] - t) < eps(FloatSpTime))
		end
	end
end

# for the shortest path from startRNode to endRNode,
# travelling with the given travel mode,
# return the index of the successor rNode of startRNode
function shortestPathNextRNode(net::Network, travelModeIndex::Int, startRNode::Int, endRNode::Int)
	assert(startRNode != endRNode)
	
	# shorthand:
	rNetTravel = net.rNetTravels[travelModeIndex]
	
	spFadjIndex = rNetTravel.spFadjIndexDelta[startRNode][endRNode]
	if spFadjIndex == 0
		spFadjIndex = net.spFadjIndexBase[startRNode, endRNode]
	end
	
	return net.rGraph.fadjList[startRNode][spFadjIndex]
end

# for the shortest path from startRNode to endRNode,
# travelling with the given travel mode,
# return the index of the rArc to travel along from startRNode
function shortestPathNextRArc(net::Network, travelModeIndex::Int, startRNode::Int, endRNode::Int)
	assert(startRNode != endRNode)
	
	# shorthand:
	rNetTravel = net.rNetTravels[travelModeIndex]
	
	spFadjIndex = rNetTravel.spFadjIndexDelta[startRNode][endRNode]
	if spFadjIndex == 0
		spFadjIndex = net.spFadjIndexBase[startRNode, endRNode]
	end
	
	return rNetTravel.spFadjArcList[startRNode][spFadjIndex]
end

# Given the travel mode, and indices of start and end node (in full graph),
# returns travel time for shortest path, and the first and last rNodes in the path (if any)
function shortestPathTravelTime(net::Network, travelModeIndex::Int, startFNode::Int, endFNode::Int)
	
	# shorthand:
	fNetTravel = net.fNetTravels[travelModeIndex]
	rNetTravel = net.rNetTravels[travelModeIndex]
	fNodeFromRNodeTime = fNetTravel.fNodeFromRNodeTime
	fNodeToRNodeTime = fNetTravel.fNodeToRNodeTime
	fNodeRArcs = net.fNodeRArcs
	fNodeToRNodes = net.fNodeToRNodes
	fNodeFromRNodes = net.fNodeFromRNodes
	isFNodeCommon = net.isFNodeCommon
	fNodeCommonFNodeIndex = net.fNodeCommonFNodeIndex
	
	# return values:
	shortestTravelTime = Inf
	rNodes = [nullIndex, nullIndex] # first and last rNodes in path
	
	if startFNode == endFNode
		shortestTravelTime = 0.0
		rNodes = [nullIndex, nullIndex]
		return shortestTravelTime, rNodes
	end
	
	# if startFNode or endFNode are common, can look up stored shortest path data
	if isFNodeCommon[startFNode]
		i = fNodeCommonFNodeIndex[startFNode]
		shortestTravelTime = fNetTravel.commonFNodeToFNodeTime[i, endFNode]
		rNodes = fNetTravel.commonFNodeToFNodeRNodes[i, endFNode]
		return shortestTravelTime, rNodes
	elseif isFNodeCommon[endFNode]
		i = fNodeCommonFNodeIndex[endFNode]
		shortestTravelTime = fNetTravel.fNodeToCommonFNodeTime[startFNode, i]
		rNodes = fNetTravel.fNodeToCommonFNodeRNodes[startFNode, i]
		return shortestTravelTime, rNodes
	end
	
	# see if travel is possible without using any nodes in rGraph
	if !isFNodeInRGraph(net, startFNode) && !isFNodeInRGraph(net, endFNode) && fNodeRArcs[startFNode] == fNodeRArcs[endFNode]
		rArc = findRArcFromFNodeToFNode(net, startFNode, endFNode)
		if rArc != nullIndex
			# calculate travel time on rArc
			fromRNode = net.rGraph.arcs[rArc].fromNodeIndex
			toRNode = net.rGraph.arcs[rArc].toNodeIndex
			travelTime = fNodeFromRNodeTime[endFNode][fromRNode] - fNodeFromRNodeTime[startFNode][fromRNode]
			assert(travelTime > 0)
			assert(abs(travelTime - (fNodeToRNodeTime[startFNode][toRNode] - fNodeToRNodeTime[endFNode][toRNode])) < eps())
			
			shortestTravelTime = travelTime # not necessarily shortest travel time, still need to check paths that use rNodes
			rNodes = [nullIndex, nullIndex]
		end
	end
	
	# for each pair of starting and ending rNodes (near given fNodes), calculate travel time
	for startRNode = fNodeToRNodes[startFNode], endRNode = fNodeFromRNodes[endFNode]
		travelTime = rNetTravel.spTimes[startRNode, endRNode]
		travelTime += fNodeToRNodeTime[startFNode][startRNode]
		travelTime += fNodeFromRNodeTime[endFNode][endRNode]
		if travelTime < shortestTravelTime
			shortestTravelTime = travelTime
			rNodes[1] = startRNode
			rNodes[2] = endRNode
		end
	end
	
	return shortestTravelTime, rNodes
end

# Find and return shortest path (as represented by list of nodes) from startFNode to endFNode,
# starting at startFNode at startTime, with a given travel mode
# Also return time at which each node in path will be reached
function shortestPath(net::Network, travelModeIndex::Int, startFNode::Int, endFNode::Int, startTime::Float)
	
	# shorthand:
	fNetTravel = net.fNetTravels[travelModeIndex]
	rNetTravel = net.rNetTravels[travelModeIndex]
	fNodeToRNodeTime = fNetTravel.fNodeToRNodeTime
	fNodeFromRNodeTime = fNetTravel.fNodeFromRNodeTime
	rArcFNodesTimes = fNetTravel.rArcFNodesTimes
	fGraph = net.fGraph
	rGraph = net.rGraph
	rArcs = rGraph.arcs
	rArcFNodes = net.rArcFNodes
	rNodeFNode = net.rNodeFNode
	fNodeRArcs = net.fNodeRArcs
	# fNodeToRNodes = net.fNodeToRNodes
	fNodeToRNodeNextFNode = net.fNodeToRNodeNextFNode
	
	# return values:
	fNodeList = Vector{Int}(0) # list of indices of nodes to travel through, in traversal order
	fNodeTimeList = Vector{Float}(0) # time at which each node in fNodeList will be reached
	
	if startFNode == endFNode
		fNodeList = [startFNode]
		fNodeTimeList = [startTime]
		return fNodeList, fNodeTimeList
	end
	
	# find which rNodes are used in the shortest path 
	(travelTime, rNodes) = shortestPathTravelTime(net, travelModeIndex, startFNode, endFNode)
	startRNode = rNodes[1]
	endRNode = rNodes[2]
	
	if startRNode == nullIndex
		assert(endRNode == nullIndex)
		# no nodes from rGraph are used, travel along a single rArc from startFNode to endFNode
		rArcIndices = fNodeRArcs[startFNode]
		assert(length(rArcIndices) <= 2)
		toRNode = rArcs[rArcIndices[1]].toNodeIndex
		if fNodeToRNodeTime[startFNode][toRNode] - fNodeToRNodeTime[endFNode][toRNode] < 0 # negative travel time
			toRNode = rArcs[rArcIndices[2]].toNodeIndex
		end
		# create fNodeList, fNodeTimeList
		fNodeList = [startFNode]
		fNodeTimeList = [startTime]
		fNode = startFNode
		currentTime = startTime
		while fNode != endFNode
			# find next arc, add next node and arc time
			nextFNode = fNodeToRNodeNextFNode[fNode][toRNode]
			arcTime = fNetTravel.arcTimes[fGraph.nodePairArcIndex[fNode, nextFNode]]
			currentTime += arcTime
			push!(fNodeList, nextFNode)
			push!(fNodeTimeList, currentTime)
			fNode = nextFNode # go to next node
		end
		return fNodeList, fNodeTimeList
	end
	
	# travel using rGraph
	if startRNode != nullIndex
		assert(endRNode != nullIndex)
		currentTime = startTime
		
		# find path from startFNode to startRNode
		fNodeList = [startFNode]
		fNodeTimeList = [startTime]
		fNode = startFNode
		startRNodeTime = startTime + fNodeToRNodeTime[startFNode][startRNode] # time to reach first rNode
		while fNode != rNodeFNode[startRNode]
			# find next arc, add next node and arc time
			nextFNode = fNodeToRNodeNextFNode[fNode][startRNode]
			# arcTime = fNodeToRNodeTime[fNode][startRNode] - fNodeToRNodeTime[nextFNode][startRNode]
			# arcTime = fNetTravel.arcTimes[fGraph.nodePairArcIndex[fNode, nextFNode]]
			# currentTime += arcTime
			currentTime = startRNodeTime - fNodeToRNodeTime[nextFNode][startRNode]
			push!(fNodeList, nextFNode)
			push!(fNodeTimeList, currentTime)
			fNode = nextFNode # go to next node
		end
		
		# find path from after startRNode to endRNode (rNodes[1] to rNodes[2])
		rNode = startRNode
		while rNode != endRNode
			# find next rArc, add fNodes
			rArc = shortestPathNextRArc(net, travelModeIndex, rNode, endRNode)
			nextRNode = rArcs[rArc].toNodeIndex
			fNodesTimesForRArc = rArcFNodesTimes[rArc]
			rNodeTime = currentTime
			fNodesOnRArc = rArcFNodes[rArc]
			i = 2
			while i <= length(fNodesOnRArc)
				fNode = fNodesOnRArc[i]
				currentTime = rNodeTime + fNodesTimesForRArc[i]
				push!(fNodeList, fNode)
				push!(fNodeTimeList, currentTime)
				i += 1 # go to next node
			end
			rNode = nextRNode # go to next node
		end
		
		# find path from after endRNode to endFNode
		fNode = rNodeFNode[endRNode]
		if fNode != endFNode
			rArcIndices = fNodeRArcs[endFNode]
			assert(length(rArcIndices) <= 2) # endFNode should only be on at most two rArcs, as it is not in rGraph
			rArc = nullIndex
			for rArcIndex in rArcIndices
				if rArcFNodes[rArcIndex][1] == fNode
					rArc = rArcIndex
				end
			end
			fNodesOnRArc = rArcFNodes[rArc]
			fNodeTime = currentTime
			i = 1
			while fNodesOnRArc[i] != endFNode
				i += 1 # go to next node on rArc
				fNode = fNodesOnRArc[i]
				currentTime = fNodeTime + fNodeFromRNodeTime[fNode][endRNode]
				push!(fNodeList, fNode)
				push!(fNodeTimeList, currentTime)
			end
		end
	end
	
	return fNodeList, fNodeTimeList
end

# set the common fNodes of the network,
# and the shortest path travel data between each fNode and commonFNodes
function setCommonFNodes!(net::Network, commonFNodes::Vector{Int})
	numFNodes = length(net.fGraph.nodes) # shorthand
	
	commonFNodes = sort(unique(commonFNodes))
	assert(length(commonFNodes) >= 1)
	assert(1 <= commonFNodes[1] && commonFNodes[end] <= numFNodes)
	numCommonFNodes = length(commonFNodes)
	
	# add common fNodes to net after calculating shortest path data
	net.isFNodeCommon = [false for i = 1:numFNodes]
	fNodeCommonFNodeIndex = [nullIndex for i = 1:numFNodes]
	for (i, commonFNode) in enumerate(commonFNodes)
		fNodeCommonFNodeIndex[commonFNode] = i
	end
	
	# calculate and store shortest path travel data between all fNodes and commonFNodes
	for fNetTravel in net.fNetTravels
		fNetTravel.commonFNodeToFNodeTime = Array{Float,2}(numCommonFNodes, numFNodes)
		fNetTravel.commonFNodeToFNodeRNodes = Array{Vector{Int},2}(numCommonFNodes, numFNodes)
		
		fNetTravel.fNodeToCommonFNodeTime = Array{Float,2}(numFNodes, numCommonFNodes)
		fNetTravel.fNodeToCommonFNodeRNodes = Array{Vector{Int},2}(numFNodes, numCommonFNodes)
		
		for commonFNode in commonFNodes, fNode = 1:numFNodes
			i = fNodeCommonFNodeIndex[commonFNode]
			
			(travelTime, rNodes) = shortestPathTravelTime(net, fNetTravel.modeIndex, commonFNode, fNode)
			fNetTravel.commonFNodeToFNodeTime[i,fNode] = travelTime
			fNetTravel.commonFNodeToFNodeRNodes[i,fNode] = rNodes
			
			(travelTime, rNodes) = shortestPathTravelTime(net, fNetTravel.modeIndex, fNode, commonFNode)
			fNetTravel.fNodeToCommonFNodeTime[fNode,i] = travelTime
			fNetTravel.fNodeToCommonFNodeRNodes[fNode,i] = rNodes
		end
	end
	
	# add common fNodes to net
	net.commonFNodes = commonFNodes
	net.fNodeCommonFNodeIndex = fNodeCommonFNodeIndex
	for commonFNode in commonFNodes
		net.isFNodeCommon[commonFNode] = true
	end
end

# Return index of the unique rArc that two fNodes have in common,
# where the rArc contains fromFNode before toFNode
# return nullIndex if no such arc found
function findRArcFromFNodeToFNode(net::Network, fromFNode::Int, toFNode::Int)
	
	assert(fromFNode != toFNode)
	assert(isFNodeInRGraph(net, fromFNode) + isFNodeInRGraph(net, toFNode) < 2) # otherwise result may not be unique, as multiple rArcs can connect two rNodes
	
	rArcIndices = intersect(net.fNodeRArcs[fromFNode], net.fNodeRArcs[toFNode])
	rArc = nullIndex
	numRArcsFound = 0
	for rArcIndex in rArcIndices
		if net.rArcFNodeIndex[rArcIndex][fromFNode] < net.rArcFNodeIndex[rArcIndex][toFNode]
			rArc = rArcIndex
			numRArcsFound += 1
		end
	end
	assert(numRArcsFound <= 1) # result should be unique
	
	return rArc
end
