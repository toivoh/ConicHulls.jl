module Hulls

export ConicHull, verify, hulltype, create_simplex_hull, add!
export ftype, gtype
export dominates
export find_dominated_facet

using ..Common
using ..Dets
import ..Common.nconic


dot(f::Facet, g::Generator) = det(g, f.generators...)
hassign(x::Number, positive::Bool) = positive ? x > 0 : x < 0
dominates(g::Generator, f::Facet) = hassign(dot(f, g), !f.positive)

function antidominates_replaced(g::Generator, f::Facet, k::Int, replacement::Generator)
    gs = copy(f.generators)
    gs[k] = replacement
    hassign(det(g, gs...), f.positive)
end


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

function find_dominated_facet(hull::ConicHull, generator::Generator)
    NF = nfacet(hull)
    facet = first(hull.facets)
    primary = first(hull.generators)
    
    while true
        if dominates(generator, facet); return facet; end
        found = false
        for k=1:NF
            if antidominates_replaced(generator, facet, k, primary)
                facet = facet.links[k]
                found = true
                break
            end
        end
        if !found; return nothing; end
    end
end

ismarked(hull::ConicHull, facet::Facet) = !(facet in hull.facets) 
mark!(hull::ConicHull, facet::Facet) = pop!(hull.facets, facet)

function mark_dominated!{F}(newfacets::Vector{(Int,F)}, hull::ConicHull{F}, 
                            generator::Generator, facet::F)
    if ismarked(hull, facet); return true; end
    if !dominates(generator, facet); return false; end

    mark!(hull, facet)
    for (k, nb) in enumerate(facet.links)
        if !mark_dominated!(newfacets, hull, generator, nb)
            # Found border: facet is dominated but not nb
            generators = copy(facet.generators)
            generators[k] = generator
            newfacet = F(generators, facet.positive)::AFacet
            newfacet.links[k] = nb
            nb.links[indexof(nb.links, facet)] = newfacet
            push!(newfacets, (k,newfacet))
        end
    end
    return true
end

function create_links!(hull::ConicHull, newfacet::Facet, k_external::Int)
    facet0 = newfacet.links[k_external]
    g0 = newfacet.generators[k_external]
    for (k, gen) in enumerate(newfacet.generators)
        if k == k_external; continue; end
        if !(newfacet.links[k] === nothing); continue; end

        lastfacet = newfacet
        gfrom = g0
        gto = newfacet.generators[k]
        facet = facet0
        while !ismarked(hull, facet)
            gnew = facet.generators[indexof(facet.links, lastfacet)]
            gto, gfrom = gnew, gto
            fnew = facet.links[indexof(facet.generators, gfrom)]
            facet, lastfacet = fnew, facet
        end
        
        newfacet.links[k] = facet
        facet.links[indexof(facet.generators, gto)] = newfacet
    end    
end

function add!{F}(hull::ConicHull{F}, generator::Generator)
    dominated = find_dominated_facet(hull, generator)
    if dominated === nothing; return false; end

    push!(hull.generators, generator)

    newfacets = Array((Int,F),0)
    mark_dominated!(newfacets, hull, generator, dominated)
    
    for (k,facet) in newfacets; create_links!(hull, facet, k); end
    for (k,facet) in newfacets; push!(hull.facets, facet); end

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
