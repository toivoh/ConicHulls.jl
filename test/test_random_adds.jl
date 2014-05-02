module TestRandomAdditions

using ConicHulls, ConicHulls.RefHull

function test(H, ngen::Int, r::Int)
    eval(ConicHulls.Hulls,:(numfacets=0))
    eval(ConicHulls.Primitives,:(numgenerators=0))
    NC, F, G = nconic(H), ftype(H), gtype(H)
        
    hull = create_simplex_hull(H)
    gs = copy(hull.generators)

    for k=1:ngen
        verify(hull)
        verify_hull(hull, gs)

        x = rand(-r:r, NC-1)
        push!(x, 1-sum(x))
        g = G(x)
        
        add!(hull, g)
        push!(gs, g)
    end
end

srand(673826715664)
for NC in 3:5
    H = hulltype(NC)
    @time for r=1:8, k=1:(1<<9); test(H, 4, r); end
end

end # module
