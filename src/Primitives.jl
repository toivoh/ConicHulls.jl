module Primitives

export AGen, AFacet

using ..Common
import ..Common.nconic


global numgenerators = 0


immutable AGen{T} <: Generator
    x::Vector{T}
    id::Int

    function AGen(x)
        global numgenerators
        g = new(T[x...], numgenerators += 1)
        g
    end
end

nconic(g::AGen) = length(g.x)

Base.getindex(g::AGen, k)  = g.x[k]
Base.show(io::IO, g::AGen) = print(io, "g$(g.id)")

end # module
