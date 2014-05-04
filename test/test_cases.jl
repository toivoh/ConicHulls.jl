module TestCases

using ConicHulls, ConicHulls.RefHull

function test(generators::Vector) 
    eval(ConicHulls.Hulls,:(numfacets=0))
    eval(ConicHulls.Primitives,:(numgenerators=0))
    test(create_simplex_hull(length(generators[1])), generators)
end
function test(hull::ConicHull, generators::Vector)
    NC, F, G = nconic(hull), ftype(hull), gtype(hull)
        
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

# Starting from a facet with the primary generator in its plane in find_dominated_facet
for k=1:100; test(Vector[ [2,-1,0], [-1,2,0] ]); end

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



end # module
