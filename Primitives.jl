module Primitives

export AGen, AFacet

using ..Common
import ..Common.nconic


global numgenerators = 0


immutable AGen{NC,T} <: Generator{NC}
    x::Vector{T}
    id::Int

    function AGen(x)
        global numgenerators
        g = new(T[x...], numgenerators += 1)
        @assert length(g.x) == NC
        g
    end
end

nconic{NC,G}(::Type{AGen{NC,G}}) = NC

Base.getindex(g::AGen, k)  = g.x[k]
Base.show(io::IO, g::AGen) = print(io, "g$(g.id)")

end # module
