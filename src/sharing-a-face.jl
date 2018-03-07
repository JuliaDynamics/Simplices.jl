function SharedFaceVolume(simplex2::Array{Float64, 2},
                        convexexp1in2::Array{Float64, 2},
                        ordered_vertices1::Vector{Int},
                        ordered_vertices2::Vector{Int})
    IntVol = 0
    n = size(simplex2, 1)

    indexextra1 = ordered_vertices1[n + 1][1]
    indexextra2 = ordered_vertices2[n + 1][1]
    convex_exp_lastvertof1 = convexexp1in2[:, indexextra1]

    extra_convex_coeff = convex_exp_lastvertof1[indexextra2]

    shared_vert_indicesin2 = ordered_vertices2[1:n]

    if extra_convex_coeff > 0
        negativecoeff = heaviside(-convex_exp_lastvertof1[shared_vert_indicesin2])
        negativeindices = round.(Int64, (negativecoeff .* shared_vert_indicesin2))
        nonnegativeindices = shared_vert_indicesin2[find(shared_vert_indicesin2 - negativeindices)]
        sigma = sum(convex_exp_lastvertof1[nonnegativeindices])
        nonnegativecoeffs = convex_exp_lastvertof1[vcat(nonnegativeindices, indexextra2)] /
                            (convex_exp_lastvertof1[indexextra2][1] + sigma)
        newpoint = simplex2[:, vcat(nonnegativeindices, indexextra2)] * nonnegativecoeffs
        intersecting_vertices = [simplex2[:, vec(shared_vert_indicesin2)] newpoint]
        IntVol = abs(det([ones(1, n + 1); intersecting_vertices]))
    end
    return IntVol
end

function SharedFaceVertices(simplex2::Array{Float64, 2},
                            convexexp1in2::Array{Float64, 2},
                            ordered_vertices1,
                            ordered_vertices2)
    n = size(simplex2, 1)
    indexextra1 = ordered_vertices1[n + 1][1]
    indexextra2 = ordered_vertices2[n + 1][1]
    convex_exp_lastvertof1 = convexexp1in2[:, indexextra1].'
    extra_convex_coeff = convex_exp_lastvertof1[indexextra2]
    shared_vert_indicesin2 = ordered_vertices2[1:n].'

    if extra_convex_coeff > 0
        negativecoeff = heaviside(-convex_exp_lastvertof1[1, shared_vert_indicesin2])
        negativeindices = round.(Int64, (negativecoeff .* shared_vert_indicesin2))
        nonnegativeindices = find(shared_vert_indicesin2 - negativeindices).'
        nonnegativeindices = shared_vert_indicesin2[nonnegativeindices]
        sigma = sum(convex_exp_lastvertof1[1, nonnegativeindices])
        nonnegativecoeffs = convex_exp_lastvertof1[1, vec([nonnegativeindices.' indexextra2])] /
                            (convex_exp_lastvertof1[indexextra2][1] + sigma)
        newpoint = simplex2[:, vec([nonnegativeindices.' [indexextra2]])] * nonnegativecoeffs
        intersecting_vertices = [simplex2[:, vec(shared_vert_indicesin2)] newpoint]
        if abs(det([ones(1, n + 1); intersecting_vertices])) > 0
            return intersecting_vertices
        else
            return Float64[]
        end
    end
end
