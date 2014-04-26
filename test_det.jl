module TestDet

import ConicHulls.Dets

for n=1:5
    for k=1:100
        A = rand(-8:8, (n,n))
        d1 = iround(det(A))
        d2 = Dets.det([A[:,k] for k=1:n]...)
        @assert d1 == d2
    end
end

end # module
