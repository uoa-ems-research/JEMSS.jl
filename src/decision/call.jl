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

# add call to list that is sorted by call priorities and then by arrival times within priorities
function addCallToQueueSortPriorityThenTime!(queuedCallList::Vector{Call}, call::Call)
	# find where to insert new call into list
	# maintain sorting by priority (lower value means a higher priority),
	# and then sorting by arrival time within priority
	# calls nearer to end of list are higher priority and arrived sooner
	i = length(queuedCallList) + 1
	while i > 1 && (call.priority > queuedCallList[i-1].priority ||
		(call.priority == queuedCallList[i-1].priority && call.arrivalTime > queuedCallList[i-1].arrivalTime))
		i -= 1
	end
	insert!(queuedCallList, i, call)
	
	if checkMode
		# check that calls are ordered first by priority, then by time
		for i = length(queuedCallList): -1 : 2
			@assert(queuedCallList[i].priority < queuedCallList[i-1].priority ||
				(queuedCallList[i].priority == queuedCallList[i-1].priority &&
				queuedCallList[i].arrivalTime <= queuedCallList[i-1].arrivalTime))
		end
	end
end

