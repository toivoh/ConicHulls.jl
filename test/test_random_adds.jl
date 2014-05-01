module TestRandomAdditions

using ConicHulls, ConicHulls.RefHull

function test(H, ngen::Int, r::Int)
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

srand(67382671)
for NC in 3:5
    H = hulltype(NC)
    @time for k=1:128; test(H, 4, 64); end
end

end # module
