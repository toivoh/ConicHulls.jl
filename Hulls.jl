module Hulls

export ConicHull

using ..Common
using ..Dets
import ..Common.nconic


type ConicHull{F,G}
    generators::Vector{G}
    facets::Set{F}
    function ConicHull(generators, facets)
        @assert nconic(F) == nconic(G)
        new(G[generators...], Set{F}(facets))
    end
end


ftype{F,G}(::Type{ConicHull{F,G}}) = F
ftype{F,G}(::ConicHull{F,G})       = F

gtype{F,G}(::Type{ConicHull{F,G}}) = G
gtype{F,G}(::ConicHull{F,G})       = G

nconic{F,G}(::Type{ConicHull{F,G}}) = nconic(G)
nconic{F,G}(::ConicHull{F,G})       = nconic(G)


function create_simplex(NC, F, G)
    gs = [G(1:NC .== k) for k in 1:NC]
    @assert det(gs...) > 0

    fs = F[F(except_index(gs, kf), isodd(kf)) for kf in 1:NC]

    for (kf, f) in enumerate(fs)
        links = except_index([1:NC...], kf)
        for (kn, nb) in enumerate(links)
            f.links[kn] = fs[nb]
        end
    end

    gs, fs
end

function create_simplex_hull{F,G}(H::Type{ConicHull{F,G}})
    NC = nconic(H)
    gs, fs = create_simplex(NC, F, G)
    @show gs
    @show fs
    ConicHull{F,G}(gs, fs)
end

end # module
