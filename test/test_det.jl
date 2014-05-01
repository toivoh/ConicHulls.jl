module TestDet

import ConicHulls.Dets

for n=1:5
    for k=1:100
        A = rand(-8:8, (n,n))
        d0 = iround(det(A))
        d1 = Dets.det([A[:,k] for k=1:n]...)
        d2 = dot(Dets.incomplete_det([A[:,k] for k=1:(n-1)]...), A[:,end])

        @assert d1 == d0
        @assert d2 == d0
    end
end

end # module
