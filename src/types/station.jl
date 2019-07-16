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

function updateStationStats!(station::Station; numIdleAmbsChange::Int = 0, time::Float = nullTime)
	if numIdleAmbsChange != 0
		@assert(time != nullTime)
		@assert(time >= station.currentNumIdleAmbsSetTime)
		
		# stats
		station.numIdleAmbsTotalDuration[station.currentNumIdleAmbs] += time - station.currentNumIdleAmbsSetTime
		station.currentNumIdleAmbs += numIdleAmbsChange
		station.currentNumIdleAmbsSetTime = time
	end
end
