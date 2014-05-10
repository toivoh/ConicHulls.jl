module TestBasic

using ConicHulls, ConicHulls.Common, ConicHulls.Dets, ConicHulls.Hulls
using ConicHulls.RefHull

function printhull(hull::ConicHull)
    println("ConicHull:")
    for facet in hull.facets
        println("  ", facet)
    end
end

function test(NC)
    initial_gs = eye(Int, NC)
    hull = create_hull(NC, initial_gs)
    NC, F, G = nconic(hull), ftype(hull), gtype(hull)
        
    gs, fs = Hulls.create_simplex(NC, F, G, initial_gs)
    # @show gs
    # @show fs
        
    @assert det(gs...) > 0
    @assert det(gs[[2, 1, (3:NC)...]]...) == -det(gs...)
    for k=1:NC
        @assert !dominates(gs[1], fs[1])
    end
        
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

    x = [-1, fill(1, NC-1)...]
#    printhull(hull)
#    println(g)
    g = add!(hull, x)
    @assert g != nothing
    verify(hull)

    push!(gs, g)
    verify_hull(hull, gs)
#    printhull(hull)
end

for k=1:2
    for NC in 3:5
        @time test(NC)
    end
end

end # module
