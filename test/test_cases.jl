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

# Errors in coplanar facter filtering
test(Vector[[1,1,0,-1]])

# Test elimination of identical (proportional) generator
test(Vector[[1,0]])
test(Vector[[2,0]])
test(Vector[[1,0,0]])
test(Vector[[2,0,0]])
test(Vector[[1,0,0,0]])
test(Vector[[2,0,0,0]])

# Test elimination of colinear generator
test(Vector[[2,-1]])
test(Vector[[2,0,-1]])
test(Vector[[2,0,0,-1]])

# Test elimination of coplanar generator
test(Vector[[-1,-1,2]])
test(Vector[[-1,-1,0,2]])

# Test four coplanar and conically independent generators
test(Vector[[1,1,-1]])
test(Vector[[1,1,-1,0]])


# Unknown issue: sometimes fail:
test(Vector[ [2,0,-1], [-1,0,2] ])
test(Vector[ [2,0,-1], [-3,0,4] ])
test(Vector[ [3,-2,0], [-1,2,0] ])

end # module
