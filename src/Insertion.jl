module Insertion

export add!

import ..Common
using ..Common
using ..Verify
using ..Hulls: ConicHull
using ..Hulls: antidominates_replaced


function find_dominated_facet(hull::ConicHull, generator::Generator)
    NF = nfacet(hull)

    # Make sure that the primary is not in the plane of the facet;
    # if the generator also is it will be unable to escape.
    # todo: more efficient way to find a primary generator and starting facet
    # that span a volume
    facet = first(hull.facets)
    primary = generatorof(facet, 1) # get a primary that is extreme
    for f in hull.facets
        if Common.dot(f, primary) != 0
            facet = f
            break
        end
    end
    @assert Common.dot(facet, primary) != 0

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

function mark_dominated!{F}(newfacets::Vector{F}, removedfacets::Vector{F}, hull::ConicHull{F}, 
                            generator::Generator, facet::F)
    if ismarked(hull, facet); return true; end
    if antidominates(generator, facet); return false; end

    mark_facet!(hull, facet)
    push!(removedfacets, facet)
    for (face, nb) in zip(facesof(facet), nbsof(facet))
        if !mark_dominated!(newfacets, removedfacets, hull, generator, nb)
            # Found border: facet is (weakly) dominated but not nb
            newfacet = add_facet!(hull, face, generator)
            mark_facet!(hull, newfacet)
            set_nb!(newfacet, nb, generator)
            set_nb!(nb, newfacet, facet)
            push!(newfacets, newfacet)
        end
    end
    return true
end

function create_links!(hull::ConicHull, newfacet::Facet, generator::Generator)
    facet0 = nbof(newfacet, generator)
    for (face, gen) in zip(facesof(newfacet), generatorsof(newfacet))
        if gen == generator; continue; end
        if !(nbof(face) === nothing); continue; end

        lastfacet = newfacet
        gfrom = generator
        gto = generatorof(face)
        facet = facet0
        while !ismarked(hull, facet)
            gnew = generatorof(facet, lastfacet)
            gto, gfrom = gnew, gto
            fnew = nbof(facet, gfrom)
            facet, lastfacet = fnew, facet
        end
        
        set_nb!(face, facet)
        set_nb!(facet, newfacet, gto)
    end    
end

add!{F,G}(hull::ConicHull{F,G}, x) = add!(hull, G(x))
function add!{F,G}(hull::ConicHull{F,G}, generator::G)
    @assert nconic(generator) == nconic(hull)
    dominated = find_dominated_facet(hull, generator)
    if dominated === nothing; return nothing; end

    push!(hull.generators, generator)

    newfacets, removedfacets = F[], F[]
    mark_dominated!(newfacets, removedfacets, hull, generator, dominated)
    
    for facet in removedfacets; del_facet!(hull, facet); end
    for facet in newfacets; create_links!(hull, facet, generator); end
    for facet in newfacets; mark_facet!(hull, facet, false); end

    return generator
end


end # module
