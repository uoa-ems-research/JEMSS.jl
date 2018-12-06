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

# EMS simulation animation

# global variables; hacky...
global animConnections = Dict{Int,WebSocket}() # store open connections
global animSimQueue = Vector{Union{Simulation,String}}() # to store sims and sim filenames between animation request and start
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

# set icons for ambulances, hospitals, etc.
function animSetIcons(client::WebSocket)
	messageDict = createMessageDict("set_icons")
	pngFileUrl(filename) = string("data:image/png;base64,", filename |> read |> base64encode)
	iconPath = joinpath(@__DIR__, "..", "..", "assets", "animation", "icons")
	icons = JSON.parsefile(joinpath(iconPath, "icons.json"))
	# set iconUrl for each icon
	for (name, icon) in icons
		icon["options"]["iconUrl"] = pngFileUrl(joinpath(iconPath, string(name, ".png")))
	end
	merge!(messageDict, icons)
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
		ambLocation = getRouteCurrentLocation!(sim.net, amb.route, sim.time)
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
	
	# determine which calls to remove, update, and add
	# need to do this after finding new ambulance locations
	# shorthand variable names:
	previousCalls = sim.previousCalls
	currentCalls = sim.currentCalls
	removeCalls = setdiff(previousCalls, currentCalls)
	updateCalls = intersect(previousCalls, currentCalls)
	addCalls = setdiff(currentCalls, previousCalls)
	changeMessageDict!(messageDict, "remove_call")
	for call in removeCalls
		call.currentLoc = Location()
		messageDict["call"] = call
		write(client, json(messageDict))
	end
	changeMessageDict!(messageDict, "move_call")
	for call in updateCalls
		updateCallLocation!(sim, call)
		messageDict["call"] = call
		write(client, json(messageDict))
	end
	changeMessageDict!(messageDict, "add_call")
	for call in addCalls
		call.currentLoc = deepcopy(call.location)
		call.movedLoc = false
		updateCallLocation!(sim, call)
		messageDict["call"] = call
		write(client, json(messageDict))
	end
	sim.previousCalls = copy(sim.currentCalls) # update previousCalls
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
	global animConnections, animSimQueue
	
	animConnections[client.id] = client
	println("Client ", client.id, " connected")
	
	# get first item in animSimQueue, or get filename now
	nextAnimItem = (length(animSimQueue) > 0 ? shift!(animSimQueue) : selectXmlFile())
	sim = nothing # init
	if typeof(nextAnimItem) == Simulation
		sim = nextAnimItem
	else
		@assert(typeof(nextAnimItem) == String)
		configFilename = nextAnimItem
		println("Initialising simulation from config:", configFilename)
		sim = initSim(configFilename; allowResim = true)
		println("...initialised")
	end
	
	# check if sim can be animated
	@assert(sim.initialised, "simulation has not been initialised; see initSim function")
	@assert(!sim.animating, "simulation is already being animated")
	@assert(!sim.writeOutput, "cannot animate simulation which is writing to output files")
	if !isdefined(sim, :backup)
		if !sim.used
			backupSim!(sim)
		else
			warn("will not be able to restart animation because the simulation has been partially run and has no backup")
		end
	end
	if sim.complete
		resetSim!(sim)
	end
	
	# set map
	messageDict = createMessageDict("set_map_view")
	messageDict["map"] = sim.map
	write(client, json(messageDict))
	
	# set sim start time
	messageDict = createMessageDict("set_start_time")
	messageDict["time"] = sim.startTime
	write(client, json(messageDict))
	
	animSetIcons(client) # set icons before adding items to map
	animAddNodes(client, sim.net.fGraph.nodes)
	animAddArcs(client, sim.net) # add first, should be underneath other objects
	animSetArcSpeeds(client, sim.map, sim.net)
	animAddBuildings(client, sim)
	animAddAmbs!(client, sim)
	
	if sim.time > sim.startTime
		# set animation to current sim state
		@assert(!sim.complete) # not sure what would happen otherwise
		sim.previousCalls = Set{Call}() # in case animation has been re-opened from same sim state
		messageDict = createMessageDict("jump_to_time")
		messageDict["time"] = sim.time
		write(client, json(messageDict))
	end
	
	sim.animating = true
	
	messageDict = createMessageDict("")
	while true
		msg = read(client) # waits for message from client
		msgString = decodeMessage(msg)
		(msgType, msgData) = parseMessage(msgString)
		
		if msgType == "prepare_next_frame"
			simTime = Float(msgData[1])
			simulateToTime!(sim, simTime)
			sim.time = simTime # otherwise sim.time only stores time of last event
			messageDict["time"] = simTime
			writeClient!(client, messageDict, "prepared_next_frame")
			
		elseif msgType == "get_next_frame"
			simTime = Float(msgData[1])
			updateFrame!(client, sim, simTime) # show updated amb locations, etc
			if !sim.complete
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
			
		elseif msgType == "update_icons"
			try
				animSetIcons(client)
			catch e
				warn("Could not update animation icons")
				warn(e)
			end
			
		elseif msgType == "disconnect"
			sim.animating = false
			close(client)
			println("Client ", client.id, " disconnected")
			break
		else
			sim.animating = false
			error("Unrecognised message: ", msgString)
		end
	end
end

"""
	function animate!(sim::Union{Simulation,Compat.Nothing} = nothing;
		configFilename::String = "", port::Int = 8001, openWindow::Bool = true)
Open a web browser window to animate the simulation.
Will animate for either `sim` or `configFilename`. If neither of these are given then there will be a prompt for the simulation configuration filename once the browser window has opened.

# Keyword arguments
- `configFilename` is the name of the configuration file to load the simulation from; can be used instead of `sim`
- `port` is the port number for the local host url, e.g. `port = 8001` will use localhost:8001; this can only be set once for all animation windows
- `openWindow` can be set to `false` to prevent the window from being opened automatically, which is useful if you wish to use a non-default browser
"""
function animate!(sim::Union{Simulation,Compat.Nothing} = nothing;
	configFilename::String = "", port::Int = 8001, openWindow::Bool = true)
	@assert(sim == nothing || configFilename == "", "can only set one of: sim, configFilename")
	global animSimQueue
	if runAnimServer(port)
		if sim != nothing
			push!(animSimQueue, sim)
		elseif configFilename != ""
			push!(animSimQueue, configFilename)
		end
		openWindow ? openLocalhost(port) : println("waiting for window with port $port to be opened")
	end
end

function animate(; configFilename::String = "", port::Int = 8001, openWindow::Bool = true)
	animate!(configFilename = configFilename, port = port, openWindow = openWindow)
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
		socket = Compat.Sockets.connect(port)
		if socket.status == 3 # = open
			println("port $port is already in use, try another")
			return false
		end
	catch
		return false
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
	if Compat.Sys.iswindows()
		run(`$(ENV["COMSPEC"]) /c start $url`)
	elseif Compat.Sys.isapple()
		run(`open $url`)
	elseif Compat.Sys.islinux() || Compat.Sys.isbsd()
		run(`xdg-open $url`)
	end
end

# opens browser window for localhost:port
function openLocalhost(port::Int)
	openUrl("http://localhost:$(port)")
end

# JSON.lower for various types, to reduce length of string returned from json function
JSON.lower(a::Ambulance) = Dict("index" => a.index, "currentLoc" => a.currentLoc)
JSON.lower(c::Call) = Dict("index" => c.index, "currentLoc" => c.currentLoc, "priority" => c.priority)
JSON.lower(h::Hospital) = Dict("index" => h.index, "location" => h.location)
JSON.lower(s::Station) = Dict("index" => s.index, "location" => s.location)
