module TestCases

using ConicHulls, ConicHulls.RefHull

test(generators::Vector) = test(hulltype(length(generators[1])), generators)
function test(H, generators::Vector)
    eval(ConicHulls.Hulls,:(numfacets=0))
    eval(ConicHulls.Primitives,:(numgenerators=0))
    NC, F, G = nconic(H), ftype(H), gtype(H)
        
    hull = create_simplex_hull(H)
    gs = copy(hull.generators)

    verify(hull)
    verify_hull(hull, gs)
    for x in generators
        g = G(x)
        add!(hull, g)
        push!(gs, g)

        verify(hull)
        verify_hull(hull, gs)
    end
end

# Starting from a facet containing the primary generator in find_dominated_facet
test(Vector[[3,0,-4]])

# Mapping g.x => g using a Dict instead of an ObjectIdDict in verify_hull
test(Vector[[1,0,0]])

end # module
