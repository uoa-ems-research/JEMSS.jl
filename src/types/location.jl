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

# copy from src to dst
function Base.copy!(dst::Location, src::Location)
	dst.x = src.x
	dst.y = src.y
end

# return true if locations have same x and y values, false otherwise
function Base.:(==)(loc1::Location, loc2::Location)
	return (loc1.x == loc2.x && loc1.y == loc2.y)
end

# return squared distance between two locations
function squareDist(map::Map, loc1::Location, loc2::Location)
	return ((loc1.x - loc2.x)*map.xScale)^2 + ((loc1.y - loc2.y)*map.yScale)^2
end

# return distance between two locations
function normDist(map::Map, loc1::Location, loc2::Location)
	return sqrt(squareDist(map, loc1, loc2))
end

# calculate off road travel time between two locations for a given travel mode
function offRoadTravelTime(travelMode::TravelMode, map::Map, loc1::Location, loc2::Location)
	return normDist(map, loc1, loc2) / travelMode.offRoadSpeed
end

# calculate off road travel time for given distance
function offRoadTravelTime(travelMode::TravelMode, dist::Float)
	return dist / travelMode.offRoadSpeed
end

# # calculate off road travel time between two locations for a given speed
# function offRoadTravelTime(map::Map, speed::Float, loc1::Location, loc2::Location)
	# return normDist(map, loc1, loc2) / speed
# end

# # calculate off road travel time for a given distance and speed
# function offRoadTravelTime(speed::Float, dist::Float)
	# return dist / speed
# end

# function linearInterp(x1::Float, x::Float, x2::Float, y1::Float, y2::Float)
	# # given x1 <= x <= x2, return matching linear interpolation for y values
	# # @assert(x1 <= x && x <= x2 && x1 < x2) # have removed this line, for speed reasons
	# return (x - x1) / (x2 - x1) * (y2 - y1) + y1
# end

function linearInterpLocation(startLoc::Location, endLoc::Location, startTime::Float, endTime::Float, currentTime::Float)
	@assert(startTime <= currentTime && currentTime <= endTime && startTime < endTime)
	# will assume constant travel time
	p = (currentTime - startTime) / (endTime - startTime) # portion of travel complete
	currentLoc = Location()
	currentLoc.x = (1-p)*startLoc.x + p*endLoc.x
	currentLoc.y = (1-p)*startLoc.y + p*endLoc.y
	
	return currentLoc
end

# return random location from map
# use uniform distribution
# allow fractional trimming of border
function randLocation(map::Map; trim::Float = 0.0, rng::AbstractRNG = GLOBAL_RNG)
	@assert(trim >= 0 && trim <= 1)
	r = rand(rng, 2) .* (1-trim) .+ trim/2
	location = Location()
	location.x = map.xMin + map.xRange * r[1]
	location.y = map.yMin + map.yRange * r[2]
	return location
end

# return map with borders trimmed by fraction
function trimmedMap(map::Map, trim::Float = 0.0)
	@assert(trim >= 0 && trim <= 1)
	mapTrimmed = deepcopy(map)
	mapTrimmed.xMin = map.xMin + map.xRange * trim / 2
	mapTrimmed.xMax = map.xMax - map.xRange * trim / 2
	mapTrimmed.yMin = map.yMin + map.yRange * trim / 2
	mapTrimmed.yMax = map.yMax - map.yRange * trim / 2
	mapTrimmed.xRange = mapTrimmed.xMax - mapTrimmed.xMin
	mapTrimmed.yRange = mapTrimmed.yMax - mapTrimmed.yMin
	return mapTrimmed
end
