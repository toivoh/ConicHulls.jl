include("ConicHulls.jl")

module Test

using ConicHulls, ConicHulls.Common, ConicHulls.Dets, ConicHulls.Hulls

for NC in 3:5

    H = hulltype(NC)
    F, G = ftype(H), gtype(H)

    gs, fs = Hulls.create_simplex(NC, F, G)
    # @show gs
    # @show fs

    @assert det(gs...) > 0
    @assert det(gs[[2, 1, (3:NC)...]]...) == -det(gs...)
    for k=1:NC
        @assert !dominates(gs[1], fs[1])
    end

    hull = create_simplex_hull(H)
    verify(hull)

    for l=1:NC
        for k=1:l
            x = gs[k].x + gs[l].x
            # Check that conic combination of generators is in the hull
            g_in = G(x) # todo: better / safer way to create new generator?
            @assert find_dominated_facet(hull, g_in) === nothing
            # Check that conic combination of surface ray and [-1,-1,...]
            # is not in the hull (since the latter isn't)
            x_out = NC*x .- 1
            g_out = G(NC*x .- 1)
            # @show x_out
            @assert !(find_dominated_facet(hull, g_out) === nothing)
        end
    end
end

end # module
