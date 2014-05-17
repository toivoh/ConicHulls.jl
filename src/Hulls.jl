module Hulls

export ConicHull, verify, validate, create_hull, init_hull!, ftype

using ..Common
using ..Dets
using ..Verify
import ..Common: nconic, gtype, get_canonical_winding, indexof
import ..Common: dot, set_opposite!, replace_link!, opposite


# --------------------------------- AFacet -----------------------------------

global numfacets = 0

immutable AFacet{G} <: Facet
    generators::Vector{G}
    links::Vector{Union(AFacet{G},Nothing)}
    positive::Bool
    id::Int

    function AFacet(generators, positive)
        global numfacets
        NF = length(generators)
        f = new(G[generators...], fill!(Array(Union(AFacet{G},Nothing), NF), nothing),
                                        positive, numfacets += 1)
        f
    end    
end
function AFacet{F<:AFacet}(face::Face{F}, generator::Generator)
    facet = face.parent
    generators = copy(facet.generators)
    generators[face.k] = generator
    F(generators, facet.positive)
end

gtype{G}( ::Type{AFacet{G}}) = G
nconic(facet::AFacet) = length(facet.generators) + 1

# Base.show(io::IO, f::AFacet) = print(io, "Facet$(f.id)")
function Base.show(io::IO, f::AFacet)
    print(io, "Facet$(f.id)(", f.generators, ", [")
    for (k, nb) in enumerate(f.links)
        print(io, "Facet$(nb.id)")
        if k < nfacet(f); print(io, ", "); end
    end
    print(io, "])");
end

indexof{F<:AFacet}(f::F, link::F)    = indexof(f.links, link)
indexof{G}(f::AFacet{G}, g::G) = indexof(f.generators, g)

opposite{F<:AFacet}(f::F, link::F) = f.generators[indexof(f, link)]
opposite{G}(f::AFacet{G}, g::G)    = f.links[     indexof(f, g)]
replace_link!{F<:AFacet}(f::F, new_link::F, link::F) = (f.links[indexof(f, link)] = new_link)
function set_opposite!{G}(f::AFacet{G}, new_link::AFacet{G}, g::G)
    f.links[indexof(f, g)] = new_link
end

function get_canonical_winding(facet::AFacet)
    perm = sortperm(facet.generators, by=g->g.id)
    generators = facet.generators[perm]
    even = facet.positive $ isoddperm(perm)
    (generators, even)
end


# -------------------------------- Dominance ---------------------------------

dot(f::AFacet, g::Generator)   = xorsign(det(f.generators..., g), !f.positive)
xorsign(x::Number, flip::Bool) = flip ? -x : x

hassign(x::Number, positive::Bool) = positive ? x > 0 : x < 0
function antidominates_replaced(g::Generator, face::Face, replacement::Generator)
    gs = copy(face.parent.generators)
    gs[face.k] = replacement
    hassign(det(gs..., g), !face.parent.positive)
end


# -------------------------------- ConicHull ---------------------------------

type ConicHull{F<:AFacet,G<:Generator}
    nc::Int
    generators::Vector{G}
    facets::Set{F}

#    ConicHull() = (@assert nconic(F) == nconic(G); new())
    ConicHull(nc::Int) = new(nc)
    function ConicHull(nc::Int, generators, facets)
        hull = ConicHull{F,G}(nc)
        init_hull!(hull, generators, facets)
        hull
    end
end

function init_hull!{F,G}(hull::ConicHull{F,G}, generators, facets)
    hull.generators = G[generators...]
    hull.facets = Set{F}(facets)
end

hulltype{G<:Generator}(::Type{G}) = ConicHull{AFacet{G},G}

ftype{F,G}(::Type{ConicHull{F,G}}) = F
ftype{F,G}(::ConicHull{F,G})       = F

gtype{F,G}(::Type{ConicHull{F,G}}) = G
gtype{F,G}(::ConicHull{F,G})       = G

nconic(hull::ConicHull) = hull.nc


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
            nb_g = opposite(nb, facet)
            # Check that the linked facets have all generators in common except the opposite ones
            link_gs = Set(facet.generators); pop!(link_gs, generator)
            link_nb_gs = Set(nb.generators); pop!(link_nb_gs, nb_g)
            @assert link_gs == link_nb_gs
            # Check that the opposite generators are different
            @assert nb_g != generator
        end
    end
end

function validate(hull::ConicHull)
    links = create_nb_dict(ftype(hull), hull.facets)
    for facet in hull.facets
        for (g,nb) in zip(facet.generators, facet.links)
            @assert nb === get_nb(links, facet, g)
        end
    end
end

create_simplex(NC, F, G, gs::Matrix) = create_simplex(NC, F, G, [gs[:,k] for k in 1:size(gs,2)])
function create_simplex(NC, F, G, gs::Vector)
    @assert length(gs) == NC
    gs = [G(x) for x in gs]
    for g in gs; (@assert nconic(g) == NC); end
    d = det(gs...)
    @assert d != 0
    odd = d < 0

    fs = F[F(except_index(gs, kf), isodd(NC-kf) $ odd) for kf in 1:NC]

    for (kf, f) in enumerate(fs)
        links = except_index([1:NC...], kf)
        for (kn, nb) in enumerate(links)
            f.links[kn] = fs[nb]
        end
    end

    gs, fs
end

create_hull{F,G}(nc::Int, H::Type{ConicHull{F,G}}) = H(nc)
function create_hull{F,G}(nc::Int, H::Type{ConicHull{F,G}}, generators)
    hull = create_hull(nc, H)
    init_hull!(hull, generators)
    hull
end

init_hull!(hull::ConicHull, gs::Matrix) = init_hull!(hull, [gs[:,k] for k in 1:size(gs,2)])
function init_hull!(hull::ConicHull, generators::Vector)
    NC = nconic(hull)
    gs, fs = create_simplex(NC, ftype(hull), gtype(hull), generators)
    init_hull!(hull, gs, fs)
    copy(gs)
end

end # module
