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

# type to contain a probability distribution sampler and corresponding random number generator
mutable struct DistrRng{T<:Sampleable}
	d::T
	rng::MersenneTwister
	
	function DistrRng(d::T, rng::MersenneTwister) where T <: Sampleable
		return new{T}(d, deepcopy(rng))
	end
	function DistrRng(d::T; seed::Int = nullIndex) where T <: Sampleable
		rng = (seed >= 0 ? MersenneTwister(seed) : MersenneTwister(rand(UInt32)))
		return new{T}(d, rng)
	end
end

global GlobalRngBackup = MersenneTwister(0); # global variable to sometimes store GLOBAL_RNG state

function copyRng!(dest::MersenneTwister, src::MersenneTwister)
	# faster than using the copy! function
	dest.seed = src.seed
	dest.state = src.state
	dest.vals = src.vals
	dest.ints = src.ints
	dest.idxF = src.idxF
	dest.idxI = src.idxI
	# # slower (but not by much):
	# for fname in [:seed, :state, :vals, :ints, :idxF, :idxI] # = fieldnames(MersenneTwister)
		# setfield!(dest, fname, getfield(src, fname))
	# end
end

# Generate next random value from distribution in distrRng, with RNG in distrRng.
# Base.rand(::AbstractRNG, ::Distribution) does not exist for all distributions, so have to use
# Base.rand(::Distribution) and after setting GLOBAL_RNG.
function Base.rand(distrRng::DistrRng, n::Int)
	rng = distrRng.rng # shorthand
	# store GLOBAL_RNG state in backup, set GLOBAL_RNG to rng
	copyRng!(GlobalRngBackup, GLOBAL_RNG)
	copyRng!(GLOBAL_RNG, rng)
	value = rand(distrRng.d, n)
	# set rng to state of GLOBAL_RNG, restore GLOBAL_RNG from backup
	copyRng!(rng, GLOBAL_RNG)
	copyRng!(GLOBAL_RNG, GlobalRngBackup) # for safety, so rng is no longer tied to GLOBAL_RNG
	return value
end
function Base.rand(distrRng::DistrRng)
	return rand(distrRng, 1)[1]
end
