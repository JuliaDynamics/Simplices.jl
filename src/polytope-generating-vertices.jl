"""
Fancy way of updating IntVert. We need to add the vertices that are common and contained.

"""

function PolytopeGeneratingVertices(Simplex1::Array{Float64, 2},
                                    Simplex2::Array{Float64, 2},
                                    IntVert::Array{Float64, 2},
                                    ConvexExpIntVert::Array{Float64, 2},
                                    convexexp1in2::Array{Float64, 2},
                                    convexexp2in1::Array{Float64, 2},
                                    ordered_vertices1::Vector{Int},
                                    ordered_vertices2::Vector{Int},
                                    numof1in2::Int,
                                    numof2in1::Int,
                                    Ncomm::Int)
    n = size(Simplex1, 1)

    numof1in2NotShared = numof1in2 - Ncomm
    numof2in1NotShared = numof2in1 - Ncomm

    # Deal first with the possibly shared vertices
    if Ncomm > 0
        #println("Deal first with the possibly shared vertices")
        IndexSharedVert1 = round.(Int64, ordered_vertices1[1:Ncomm].')
        IndexSharedVert2 = round.(Int64, ordered_vertices2[1:Ncomm].')
        IntVert = [IntVert; Simplex1[:, vec(IndexSharedVert1)].']

        common_vertices_expressed_in1 = convexexp2in1[:, vec(IndexSharedVert2)].'
        common_vertices_expressed_in2 = convexexp1in2[:, vec(IndexSharedVert1)].'

        ConvexExpIntVert = [ConvexExpIntVert; [common_vertices_expressed_in1 common_vertices_expressed_in2]]
    end

    # Dealing with vertices that are contained but not shared

    if numof1in2NotShared * numof2in1NotShared > 0
        # Each simplex contains vertices of the other beyond the shared ones
        #println("# Each simplex contains vertices of the other beyond the shared ones")
        Id = eye(n+1, n+1)

        indexvert1in2notshared = ordered_vertices1[Ncomm+1:numof1in2].'
        indexvert2in1notshared = ordered_vertices2[Ncomm+1:numof2in1].'

        IntVert = vcat(IntVert,
                    Simplex1[:, vec(indexvert1in2notshared)].',
                    Simplex2[:, vec(indexvert2in1notshared)].')

        conv_exp_vertof1_expressed_in1 = Id[vec(indexvert1in2notshared), :]
        conv_exp_vertof2_expressed_in2 = Id[vec(indexvert2in1notshared), :]
        conv_exp_vertof1_expressed_in2 = convexexp1in2[:, vec(indexvert1in2notshared)].'
        conv_exp_vertof2_expressed_in1 = convexexp2in1[:, vec(indexvert2in1notshared)].'

        ConvexExpIntVert = vcat(ConvexExpIntVert,
                                    [conv_exp_vertof1_expressed_in1 conv_exp_vertof1_expressed_in2],
                                    [conv_exp_vertof2_expressed_in1 conv_exp_vertof2_expressed_in2])

    elseif numof1in2NotShared > 0
        # Only Simplex2 contains vertices of Simplex1 other than the shared ones
        #println("Only Simplex2 contains vertices of Simplex1 other than the shared ones")
        Id = eye(n+1, n+1)

        indexvert1in2notshared = ordered_vertices1[Ncomm+1:numof1in2].'

        IntVert = vcat(IntVert,
                    Simplex1[:, vec(indexvert1in2notshared)].')

        conv_exp_vertof1_expressed_in1 = Id[vec(indexvert1in2notshared), :]
        conv_exp_vertof1_expressed_in2 = convexexp1in2[:, vec(indexvert1in2notshared)].'

        ConvexExpIntVert = vcat(ConvexExpIntVert,
                                    [conv_exp_vertof1_expressed_in1 conv_exp_vertof1_expressed_in2])
    elseif numof2in1NotShared > 0
        # Only Simplex1 contains vertices of Simplex2 other than the shared ones
        #println("Only Simplex1 contains vertices of Simplex2 other than the shared ones")

        Id = eye(n+1, n+1)

        indexvert2in1notshared = ordered_vertices2[Ncomm+1:numof2in1].'

        IntVert = vcat(IntVert,
                    Simplex2[:, vec(indexvert2in1notshared)].')

        conv_exp_vertof2_expressed_in2 = Id[vec(indexvert2in1notshared), :]
        conv_exp_vertof2_expressed_in1 = convexexp2in1[:, vec(indexvert2in1notshared)].'
        ConvexExpIntVert = vcat(ConvexExpIntVert,
                                    [conv_exp_vertof2_expressed_in1 conv_exp_vertof2_expressed_in2])
    end

    return IntVert, ConvexExpIntVert
end
