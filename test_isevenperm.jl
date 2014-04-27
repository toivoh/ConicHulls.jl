module TestIsevenperm

using ConicHulls.Common.isevenperm, ConicHulls.Dets.det

for n=1:5
    eyerows = [int(1:n .== k) for k=1:n]
    for p in permutations(1:n)
        even = isevenperm(p)
        d = det(eyerows[p]...)
        @assert even ? d == 1 : d == -1
    end
end

end # module
