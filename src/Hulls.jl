module Hulls

export ConicHull, verify, validate, create_hull, init_hull!, add!, ftype

using ..Common
using ..Dets
using ..Verify
import ..Common: nconic, gtype, get_canonical_winding, indexof


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

Base.length(faces::Faces) = nfacet(faces.parent)
Base.start(faces::Faces) = 1
Base.next{F}(faces::Faces{F}, k) = (Face{F}(faces.parent, k), k+1)
Base.done(faces::Faces, k) = k == (nfacet(faces.parent)+1)
Base.eltype{F}(faces::Faces{F}) = Face{F}

facesof(f::Facet) = Faces(f)


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

dot(f::Facet, g::Generator) = det(f.generators..., g)
hassign(x::Number, positive::Bool) = positive ? x > 0 : x < 0
dominates(g::Generator, f::Facet) = hassign(dot(f, g), f.positive)
antidominates(g::Generator, f::Facet) = hassign(dot(f, g), !f.positive)

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

function find_dominated_facet(hull::ConicHull, generator::Generator)
    NF = nfacet(hull)

    # Make sure that the primary is not in the plane of the facet;
    # if the generator also is it will be unable to escape.
    # todo: more efficient way to find a primary generator and starting facet
    # that span a volume
    facet = first(hull.facets)
    primary = facet.generators[1] # get a primary that is extreme
    for f in hull.facets
        if dot(f, primary) != 0
            facet = f
            break
        end
    end
    @assert dot(facet, primary) != 0

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
    if antidominates(generator, facet); return false; end

    mark!(hull, facet)
    for (face, nb) in zip(facesof(facet), facet.links)
        if !mark_dominated!(newfacets, hull, generator, nb)
            # Found border: facet is (weakly) dominated but not nb
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
    for (face, gen) in zip(facesof(newfacet), newfacet.generators)
        if gen == generator; continue; end
        if !(nbof(face) === nothing); continue; end

        lastfacet = newfacet
        gfrom = generator
        gto = opposite(face)
        facet = facet0
        while !ismarked(hull, facet)
            gnew = opposite(facet, lastfacet)
            gto, gfrom = gnew, gto
            fnew = opposite(facet, gfrom)
            facet, lastfacet = fnew, facet
        end
        
        set_nb!(face, facet)
        set_opposite!(facet, newfacet, gto)
    end    
end

add!{F,G}(hull::ConicHull{F,G}, x) = add!(hull, G(x))
function add!{F,G}(hull::ConicHull{F,G}, generator::G)
    @assert nconic(generator) == nconic(hull)
    dominated = find_dominated_facet(hull, generator)
    if dominated === nothing; return nothing; end

    push!(hull.generators, generator)

    newfacets = F[]
    mark_dominated!(newfacets, hull, generator, dominated)
    
    for facet in newfacets; create_links!(hull, facet, generator); end
    for facet in newfacets; push!(hull.facets, facet); end

    return generator
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
