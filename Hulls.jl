module Hulls

export ConicHull, verify, hulltype, create_simplex_hull
export ftype, gtype
export dominates

using ..Common
using ..Dets
import ..Common.nconic


dot(f::Facet, g::Generator) = det(g, f.generators...)
dominates(g::Generator, f::Facet) = (dot(f, g) < 0) == f.positive


# --------------------------------- AFacet -----------------------------------

global numfacets = 0

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


# -------------------------------- ConicHull ---------------------------------

type ConicHull{F<:AFacet,G<:Generator}
    generators::Vector{G}
    facets::Set{F}
    function ConicHull(generators, facets)
        @assert nconic(F) == nconic(G)
        new(G[generators...], Set{F}(facets))
    end
end

hulltype{G<:Generator}(::Type{G}) = ConicHull{AFacet{nfacet(G),G},G}

ftype{F,G}(::Type{ConicHull{F,G}}) = F
ftype{F,G}(::ConicHull{F,G})       = F

gtype{F,G}(::Type{ConicHull{F,G}}) = G
gtype{F,G}(::ConicHull{F,G})       = G

nconic{F,G}(::Type{ConicHull{F,G}}) = nconic(G)
nconic{F,G}(::ConicHull{F,G})       = nconic(G)


function verify(hull::ConicHull)
    NF = nfacet(hull)
    for facet in hull.facets
        @assert length(facet.links) == length(facet.generators) == NF
        if length(Set(facet.links)) != NF
            error("facet $facet has duplicate neighbors") 
        end
        if length(Set(facet.generators)) != NF
            error("facet $facet has duplicate generators") 
        end
        for (nb, generator) in zip(facet.links, facet.generators)
            # @assert !is(nb, nothing)
            if !(nb in hull.facets)
                error("facet $facet links to facet $nb, which is not in the hull.")
            end
            # Check that neighbor is not dominated by the opposing generator
            @assert !dominates(generator, nb)
            
            # Check that there is a link back
            l = indexof(nb.links, facet)
            # Check that the linked facets have all generators in common except the opposite ones
            link_gs = Set(facet.generators); pop!(link_gs, generator)
            link_nb_gs = Set(except_index(nb.generators, l));
            @assert link_gs == link_nb_gs
            # Check that the opposite generators are different
            @assert (nb.generators[l] != generator)
        end
    end
end

function create_simplex(NC, F, G)
    gs = [G(1:NC .== k) for k in 1:NC]
    @assert det(gs...) > 0

    fs = F[F(except_index(gs, kf), isodd(kf)) for kf in 1:NC]

    for (kf, f) in enumerate(fs)
        links = except_index([1:NC...], kf)
        for (kn, nb) in enumerate(links)
            f.links[kn] = fs[nb]
        end
    end

    gs, fs
end

function create_simplex_hull{F,G}(H::Type{ConicHull{F,G}})
    NC = nconic(H)
    gs, fs = create_simplex(NC, F, G)
    ConicHull{F,G}(gs, fs)
end

end # module
