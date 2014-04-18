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
end

end # module
