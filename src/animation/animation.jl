# ambulance simulation animation

# global variables; hacky...
global animConnections = Dict{Int,WebSocket}() # store open connections
global animConfigFilenames = Vector{String}() # store filenames between animation request and start
global animPort = nullIndex # localhost port for animation, to be set

function decodeMessage(msg)
	return String(msg)
end

# parse message from html, extract values
function parseMessage(msg::String)
	msgSplit = readdlm(IOBuffer(msg))
	# msgType = msgSplit[1]
	# msgData = msgSplit[2:end]
	return msgSplit[1], msgSplit[2:end]
end

# create dictionary for sending messages to js
function createMessageDict(message::String)
	messageDict = Dict()
	messageDict["message"] = message
	return messageDict
	# return Dict([("message", message)])
end

# change message of messageDict
function changeMessageDict!(messageDict::Dict, message::String)
	messageDict["message"] = message
end

function writeClient!(client::WebSocket, messageDict::Dict, message::String)
	# common enough lines to warrant the use of a function, I guess
	changeMessageDict!(messageDict, message)
	write(client, json(messageDict))
end

# adds nodes from fGraph
function animAddNodes(client::WebSocket, nodes::Vector{Node})
	messageDict = createMessageDict("add_node")
	for node in nodes
		messageDict["node"] = node
		write(client, json(messageDict))
	end
end

# adds arcs from rGraph, should only be called after animAddNodes()
function animAddArcs(client::WebSocket, net::Network)
	messageDict = createMessageDict("add_arc")
	for arc in net.rGraph.arcs
		messageDict["arc"] = arc
		messageDict["node_indices"] = net.rArcFNodes[arc.index]
		write(client, json(messageDict))
	end
end

# calculate and set speeds for arcs in rGraph, should only be called after animAddArcs()
function animSetArcSpeeds(client::WebSocket, map::Map, net::Network)
	# shorthand:
	fNodes = net.fGraph.nodes 
	rNetTravels = net.rNetTravels
	
	messageDict = createMessageDict("set_arc_speed")
	for arc in net.rGraph.arcs
		# calculate arc distance
		dist = 0.0
		fNodeIndices = net.rArcFNodes[arc.index]
		for i = 1:length(fNodeIndices)-1
			dist += normDist(map, fNodes[fNodeIndices[i]].location, fNodes[fNodeIndices[i+1]].location)
		end
		
		# select minimum travel time
		arcTime = Inf
		for i = 1:length(rNetTravels)
			if arcTime > rNetTravels[i].arcTimes[arc.index]
				arcTime = rNetTravels[i].arcTimes[arc.index]
			end
		end
		
		messageDict["arc"] = arc
		messageDict["speed"] = dist / (arcTime*24) # km / hr
		write(client, json(messageDict))
	end
end

function animAddBuildings(client::WebSocket, sim::Simulation)
	messageDict = createMessageDict("add_hospital")
	for h in sim.hospitals
		messageDict["hospital"] = h
		write(client, json(messageDict))
	end
	delete!(messageDict, "hospital")
	
	changeMessageDict!(messageDict, "add_station")
	for s in sim.stations
		messageDict["station"] = s
		write(client, json(messageDict))
	end
	# delete!(messageDict, "station")
end

function animAddAmbs!(client::WebSocket, sim::Simulation)
	messageDict = createMessageDict("add_ambulance")
	for amb in sim.ambulances
		ambLocation = getRouteCurrentLocation!(sim.net, amb.route, sim.startTime)
		amb.currentLoc = ambLocation
		messageDict["ambulance"] = amb
		write(client, json(messageDict))
	end
end

# write frame updates to client
function updateFrame!(client::WebSocket, sim::Simulation, time::Float)
	
	# check which ambulances have moved since last frame
	# need to do this before showing call locations
	messageDict = createMessageDict("move_ambulance")
	for amb in sim.ambulances
		ambLocation = getRouteCurrentLocation!(sim.net, amb.route, time)
		if !isSameLocation(ambLocation, amb.currentLoc)
			amb.currentLoc = ambLocation
			amb.movedLoc = true
			# move ambulance
			messageDict["ambulance"] = amb
			write(client, json(messageDict))
		else
			amb.movedLoc = false
		end
	end
	delete!(messageDict, "ambulance")
	
	# display current calls by comparing sim.previousCallList and sim.currentCallList
	# need to do this after finding new ambulance locations
	# shorthand variable names:
	previousCallList = sim.previousCallList
	currentCallList = sim.currentCallList
	numPreviousCalls = length(previousCallList)
	numCurrentCalls = length(currentCallList)
	if checkMode
		# code expects call lists to remain sorted by call indices
		assert(all(previousCallList[i].index < previousCallList[i+1].index for i = 1:numPreviousCalls - 1))
		assert(all(currentCallList[i].index < currentCallList[i+1].index for i = 1:numCurrentCalls - 1))
	end
	changeMessageDict!(messageDict, "")
	j = 1 # for indexing currentCallList
	for i = 1:numPreviousCalls
		if j <= numCurrentCalls && previousCallList[i].index == currentCallList[j].index
			# call still exists, consider moving it if call status indicates location other than call origin
			call = currentCallList[j]
			updateCallLocation!(sim, call)
			changeMessageDict!(messageDict, "move_call")
			messageDict["call"] = call
			write(client, json(messageDict))
			# move to comparing next call in current call list
			j += 1
		else
			# previous call does not exist in current call list, remove
			call = previousCallList[i]
			call.currentLoc = Location()
			changeMessageDict!(messageDict, "remove_call")
			messageDict["call"] = call
			write(client, json(messageDict))
		end
	end
	changeMessageDict!(messageDict, "add_call")
	for i = j:numCurrentCalls
		# new call - set current location to origin, move if needed
		call = currentCallList[i]
		call.currentLoc = deepcopy(call.location)
		call.movedLoc = false
		updateCallLocation!(sim, call)
		messageDict["call"] = call
		write(client, json(messageDict))
	end
	# update previousCallList
	sim.previousCallList = copy(sim.currentCallList)
end

# update call current location
function updateCallLocation!(sim::Simulation, call::Call)
	# consider moving call if the status indicates location other than call origin location
	if call.status == callGoingToHospital || call.status == callAtHospital
		amb = sim.ambulances[call.ambIndex]
		call.movedLoc = amb.movedLoc
		if amb.movedLoc
			call.currentLoc = amb.currentLoc
		end
	end
end

wsh = WebSocketHandler() do req::Request, client::WebSocket
	global animConnections, animConfigFilenames
	
	animConnections[client.id] = client
	println("Client ", client.id, " connected")
	
	# get oldest filename from animConfigFilenames, or select file now
	configFilename = (length(animConfigFilenames) > 0 ? shift!(animConfigFilenames) : selectXmlFile())
	println("Running from config: ", configFilename)
	
	println("Initialising simulation...")
	sim = initSimulation(configFilename; allowResim = true)
	println("...initialised")
	
	# set map
	messageDict = createMessageDict("set_map_view")
	messageDict["map"] = sim.map
	write(client, json(messageDict))
	
	# set sim start time
	messageDict = createMessageDict("set_start_time")
	messageDict["time"] = sim.startTime
	write(client, json(messageDict))
	
	animAddNodes(client, sim.net.fGraph.nodes)
	animAddArcs(client, sim.net) # add first, should be underneath other objects
	animSetArcSpeeds(client, sim.map, sim.net)
	animAddBuildings(client, sim)
	animAddAmbs!(client, sim)
	
	messageDict = createMessageDict("")
	while true
		msg = read(client) # waits for message from client
		msgString = decodeMessage(msg)
		(msgType, msgData) = parseMessage(msgString)
		
		if msgType == "prepare_next_frame"
			simTime = Float(msgData[1])
			simEnd = simulateToTime!(sim, simTime)
			messageDict["time"] = simTime
			writeClient!(client, messageDict, "prepared_next_frame")
			
		elseif msgType == "get_next_frame"
			simTime = Float(msgData[1])
			updateFrame!(client, sim, simTime) # show updated amb locations, etc
			if length(sim.eventList) > 0
				messageDict["time"] = simTime
				writeClient!(client, messageDict, "got_next_frame")
			else
				# no events left, finish animation
				writeClient!(client, messageDict, "got_last_frame")
			end
			
		elseif msgType == "pause"
		
		elseif msgType == "stop"
			# reset
			resetSim!(sim)
			animAddAmbs!(client, sim)
			
		elseif msgType == "disconnect"
			close(client)
			println("Client ", client.id, " disconnected")
			break
		else
			error("Unrecognised message: ", msgString)
		end
	end
end

# start the animation, open a browser window for it
# can set the port for the connection, and the simulation config filename
# openWindow = false prevents a browser window from opening automatically, will need to open manually
function animate(; port::Int = 8001, configFilename::String = "", openWindow::Bool = true)
	global animConfigFilenames
	if runAnimServer(port)
		if configFilename != ""
			push!(animConfigFilenames, configFilename)
		end
		openWindow ? openLocalhost(port) : println("waiting for window with port $port to be opened")
	end
end

# creates and runs server for given port
# returns true if server is running, false otherwise
function runAnimServer(port::Int)
	# check if port already in use
	global animPort
	if port == animPort && port != nullIndex
		return true # port already used for animation
	elseif animPort != nullIndex
		println("use port $animPort instead")
		return false
	end
	try
		socket = connect(port)
		if socket.status == 3 # = open
			println("port $port is already in use, try another")
			return false
		end
	end
	
	# create and run server
	onepage = readstring("$sourcePath/animation/index.html")
	httph = HttpHandler() do req::Request, res::Response
		Response(onepage)
	end
	server = Server(httph, wsh)
	@async run(server, port)
	animPort = port
	println("opened port $animPort, use this for subsequent animation windows")
	return true
end

# opens browser window for url
function openUrl(url::String)
	if is_windows()
		run(`$(ENV["COMSPEC"]) /c start $url`)
	elseif is_apple()
		run(`open $url`)
	elseif is_linux() || is_bsd()
		run(`xdg-open $url`)
	end
end

# opens browser window for localhost:port
function openLocalhost(port::Int)
	openUrl("http://localhost:$(port)")
end
