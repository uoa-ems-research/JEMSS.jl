# type to contain a probability distribution and corresponding random number generator
type DistrRng
	d::Distribution
	rng::MersenneTwister
	
	# DistrRng(d, rng) = new(r, deepcopy(rng))
	function DistrRng(d::Distribution; seed::Int = nullIndex)
		rng = (seed >= 0 ? MersenneTwister(seed) : MersenneTwister(rand(UInt32)))
		return new(d, rng)
	end
end

global GlobalRngBackup = MersenneTwister(0); # global variable to sometimes store GLOBAL_RNG state

function copyRng!(dest::MersenneTwister, src::MersenneTwister)
	# faster than using the copy! function
	dest.seed = src.seed
	dest.state = src.state
	dest.vals = src.vals
	dest.idx = src.idx
	# # slower (but not by much):
	# for fname in [:seed, :state, :vals, :idx] # = fieldnames(MersenneTwister)
		# setfield!(dest, fname, getfield(src, fname))
	# end
end

# Generate next random value from distribution in distrRng, with RNG in distrRng.
# Base.rand(::AbstractRNG, ::Distribution) does not exist for all distributions, so have to use
# Base.rand(::Distribution) and after setting GLOBAL_RNG.
function Base.rand(distrRng::DistrRng)
	rng = distrRng.rng # shorthand
	# store GLOBAL_RNG state in backup, set GLOBAL_RNG to rng
	copyRng!(GlobalRngBackup, Base.GLOBAL_RNG)
	copyRng!(Base.GLOBAL_RNG, rng)
	value = rand(distrRng.d)
	# set rng to state of GLOBAL_RNG, restore GLOBAL_RNG from backup
	copyRng!(rng, Base.GLOBAL_RNG)
	copyRng!(Base.GLOBAL_RNG, GlobalRngBackup) # for safety, so rng is no longer tied to GLOBAL_RNG
	return value
end
