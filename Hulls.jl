module Hulls

export ConicHull, verify, hulltype, create_simplex_hull, add!
export ftype, gtype
export dominates
export facesof
export find_dominated_facet

using ..Common
using ..Dets
import ..Common.nconic

import ..Common.indexof


# ---------------------------------- Face ------------------------------------

immutable Face{F<:Facet}
    parent::F
    k::Int
end
Face{F<:Facet}(parent::F, nb::F) = Face{F}(parent, indexof(parent, nb))
Face{F<:Facet}(parent::F, g::Generator) = Face{F}(parent, indexof(parent, g))

parentof(f::Face) = f.parent
nbof(f::Face)     = f.parent.links[f.k]
opposite(f::Face) = f.parent.generators[f.k]
flip(f::Face)     = (nb = nbof(f); Face(nb, indexof(nb, f)))

set_nb!{F}(f::Face{F}, newnb::F) = (f.parent.links[f.k] = newnb)


# ---------------------------------- Faces -----------------------------------

immutable Faces{F<:Facet}; parent::F; end

Base.length{F}(faces::Faces{F}) = nfacet(F)
Base.start(faces::Faces) = 1
Base.next{F}(faces::Faces{F}, k) = (Face{F}(faces.parent, k), k+1)
Base.done{F}(faces::Faces{F}, k) = k == (nfacet(F)+1)
Base.eltype{F}(faces::Faces{F}) = Face{F}

facesof(f::Facet) = Faces(f)


# --------------------------------- AFacet -----------------------------------

global numfacets = 0

immutable AFacet{NF,G} <: Facet{NF}
    generators::Vector{G}
    links::Vector{Union(AFacet{NF,G},Nothing)}
    positive::Bool
    id::Int

    function AFacet(generators, positive)
        global numfacets
        f = new(G[generators...], fill!(Array(Union(AFacet{NF,G},Nothing), NF), nothing),
                                        positive, numfacets += 1)
        @assert length(f.generators) == NF
        f
    end    
end
function AFacet{F<:AFacet}(face::Face{F}, generator::Generator)
    facet = face.parent
    generators = copy(facet.generators)
    generators[face.k] = generator
    F(generators, facet.positive)
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

indexof{F<:AFacet}(f::F, link::F)    = indexof(f.links, link)
indexof{NF,G}(f::AFacet{NF,G}, g::G) = indexof(f.generators, g)

opposite{F<:AFacet}(f::F, link::F)    = f.generators[indexof(f, link)]
opposite{NF,G}(f::AFacet{NF,G}, g::G) = f.links[     indexof(f, g)]
replace_link!{F<:AFacet}(f::F, new_link::F, link::F) = (f.links[indexof(f, link)] = new_link)
function set_opposite!{NF,G}(f::AFacet{NF,G}, new_link::AFacet{NF,G}, g::G)
    f.links[indexof(f, g)] = new_link
end


# -------------------------------- Dominance ---------------------------------

dot(f::Facet, g::Generator) = det(g, f.generators...)
hassign(x::Number, positive::Bool) = positive ? x > 0 : x < 0
dominates(g::Generator, f::Facet) = hassign(dot(f, g), !f.positive)

function antidominates_replaced(g::Generator, face::Face, replacement::Generator)
    gs = copy(face.parent.generators)
    gs[face.k] = replacement
    hassign(det(g, gs...), face.parent.positive)
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

function find_dominated_facet(hull::ConicHull, generator::Generator)
    NF = nfacet(hull)
    facet = first(hull.facets)
    primary = first(hull.generators)
    
    while true
        if dominates(generator, facet); return facet; end
        found = false
        for face in facesof(facet)
            if antidominates_replaced(generator, face, primary)
                facet = nbof(face)
                found = true
                break
            end
        end
        if !found; return nothing; end
    end
end

ismarked(hull::ConicHull, facet::Facet) = !(facet in hull.facets) 
mark!(hull::ConicHull, facet::Facet) = pop!(hull.facets, facet)

function mark_dominated!{F}(newfacets::Vector{F}, hull::ConicHull{F}, 
                            generator::Generator, facet::F)
    if ismarked(hull, facet); return true; end
    if !dominates(generator, facet); return false; end

    mark!(hull, facet)
    for (face, nb) in zip(facesof(facet), facet.links)
        if !mark_dominated!(newfacets, hull, generator, nb)
            # Found border: facet is dominated but not nb
            newfacet = AFacet(face, generator)
            set_opposite!(newfacet, nb, generator)
            replace_link!(nb, newfacet, facet)
            push!(newfacets, newfacet)
        end
    end
    return true
end

function create_links!(hull::ConicHull, newfacet::Facet, generator::Generator)
    facet0 = opposite(newfacet, generator)
    for (k, gen) in enumerate(newfacet.generators)
        if gen == generator; continue; end
        if !(newfacet.links[k] === nothing); continue; end

        lastfacet = newfacet
        gfrom = generator
        gto = newfacet.generators[k]
        facet = facet0
        while !ismarked(hull, facet)
            gnew = opposite(facet, lastfacet)
            gto, gfrom = gnew, gto
            fnew = opposite(facet, gfrom)
            facet, lastfacet = fnew, facet
        end
        
        newfacet.links[k] = facet
        set_opposite!(facet, newfacet, gto)
    end    
end

function add!{F}(hull::ConicHull{F}, generator::Generator)
    dominated = find_dominated_facet(hull, generator)
    if dominated === nothing; return false; end

    push!(hull.generators, generator)

    newfacets = F[]
    mark_dominated!(newfacets, hull, generator, dominated)
    
    for facet in newfacets; create_links!(hull, facet, generator); end
    for facet in newfacets; push!(hull.facets, facet); end

    return true
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
