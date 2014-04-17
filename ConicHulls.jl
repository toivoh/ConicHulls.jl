module ConicHulls

include("Common.jl")
include("Dets.jl")
include("Primitives.jl")
include("Hulls.jl")

using .Primitives
using .Hulls


function hulltype(NC)
    G = AGen{NC,Int}
    F = AFacet{NC-1,G}
    ConicHull{F,G}
end

function create_simplex_hull(NC::Int)
    Hulls.create_simplex_hull(hulltype(NC))
end


end # module
