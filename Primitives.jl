module Primitives

export AGen, AFacet

using ..Common
import ..Common.nconic


global numgenerators = 0
global numfacets = 0


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

immutable AFacet{NF,G} <: Facet{NF}
    generators::Vector{G}
    links::Vector{AFacet{NF,G}}
    positive::Bool
    id::Int

    function AFacet(generators, positive)
        global numfacets
        f = new(G[generators...], Array(AFacet{NF,G}, NF), positive, numfacets += 1)
        @assert length(f.generators) == NF
        f
    end    
end

nconic{NF,G}(::Type{AFacet{NF,G}}) = NF + 1

# Base.show(io::IO, f::AFacet) = print(io, "Facet$(f.id)")
function Base.show{NF,G}(io::IO, f::AFacet{NF,G})
    print(io, "Facet$(f.id)(", f.generators, ", [")
    for (k, nb) in enumerate(f.links)
        print(io, "Facet$(nb.id)")
        if k < NF; print(io, ", "); end
    end
    print(io, "])");
end


end # module
