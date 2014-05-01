include("../src/ConicHulls.jl")
include("../src/RefHull.jl")

module Test

using ConicHulls, ConicHulls.Common, ConicHulls.Dets, ConicHulls.Hulls
using RefHull

function printhull(hull::ConicHull)
    println("ConicHull:")
    for facet in hull.facets
        println("  ", facet)
    end
end

function verify_hull(hull::ConicHull, generators::Vector) 
    validate(hull)

    gmap = [g.x => g for g in generators]
    generators = [g.x for g in generators]
    generators, facets = RefHull.conichull(nconic(hull), generators)

    @assert length(hull.facets) == length(facets)
    facetset = Set()
    for facet in facets
        push!(facetset, tuple([gmap[g] for g in facet.generators]..., facet.positive))
    end
    for facet in hull.facets
        gs, even = get_canonical_winding(facet)
        key = tuple(gs..., even)
        @assert key in facetset
    end
end

function test(H)
    NC, F, G = nconic(H), ftype(H), gtype(H)
        
    gs, fs = Hulls.create_simplex(NC, F, G)
    # @show gs
    # @show fs
        
    @assert det(gs...) > 0
    @assert det(gs[[2, 1, (3:NC)...]]...) == -det(gs...)
    for k=1:NC
        @assert !dominates(gs[1], fs[1])
    end
        
    hull = create_simplex_hull(H)
    gs = copy(hull.generators)
    verify(hull)
    verify_hull(hull, gs)
        
    for facet in hull.facets
        x = sum([g.x for g in facet.generators])
        # Check that conic combination of generators is in the hull
        g_in = G(x) # todo: better / safer way to create new generator?
        @assert find_dominated_facet(hull, g_in) === nothing
        # Check that conic combination of the facet generators minus
        # the opposing generator dominates the facet (and no other)
        opposite_g = first(setdiff(Set(hull.generators), Set(facet.generators)))
        x_out = x - opposite_g.x
        g_out = G(x_out)
        # @show x_out
        @assert find_dominated_facet(hull, g_out) === facet
    end

    g = G([-1, fill(1, NC-1)...])
#    printhull(hull)
#    println(g)
    @assert add!(hull, g)
    verify(hull)

    push!(gs, g)
    verify_hull(hull, gs)
#    printhull(hull)
end

for k=1:2
    for NC in 3:5
        H = hulltype(NC)
        @time test(H)
    end
end

end # module
