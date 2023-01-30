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

"""
    flatten(dict::Dict{String,T}; delim = "_") where T <: Any
Given a dict of nested dicts with string keys, return a dict that is flattened (not nested) with keys between root and leaf concatenated by `delim`.

# Example
```julia-repl
julia> flatten(Dict("a" => 1, "b" => Dict("c" => 2)))
Dict{String,Any} with 2 entries:
  "b_c" => 2
  "a"   => 1
```
"""
function flatten(dict::Dict{String,T}; delim="_") where {T<:Any}
    dflat = Dict{String,Any}()
    function recurse(d; s=nothing)
        if isa(d, Dict{String,T} where {T<:Any})
            for (k, v) in d
                recurse(v; s=s == nothing ? k : string(s, delim, k))
            end
        else
            dflat[s] = d
        end
    end
    recurse(dict)
    return dflat
end
