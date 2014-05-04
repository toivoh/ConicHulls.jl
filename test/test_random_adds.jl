module TestRandomAdditions

using ConicHulls, ConicHulls.RefHull

function test(NC::Int, ngen::Int, r::Int)
    eval(ConicHulls.Hulls,:(numfacets=0))
    eval(ConicHulls.Primitives,:(numgenerators=0))
    test(create_hull(NC), ngen, r)
end

function test(hull::ConicHull, ngen::Int, r::Int)
    NC, F, G = nconic(hull), ftype(hull), gtype(hull)
    
    gs = init_hull!(hull, eye(Int, NC))

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
    @time for r=1:8, ngen=1:8, k=1:(1<<(8-NC))/ngen; test(NC, ngen, r); end
end

end # module
