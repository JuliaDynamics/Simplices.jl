
"""
	  simplexintersection(S1::Array{Float64, 2}, S2::Array{Float64, 2};
        tolerance::Float64 = 1/10^10) -> Float64

Computes the volume of intersection between two n-dimensional simplices
by boundary triangulation. The simplices `S1` and `S2` are arrays of
(n, n+1), where each column is a vertex.

## How are intersections computed?
Intersections are computed as follows:

1. Find minimal set of points generating the intersection volume. These points form
a convex polytope Pᵢ.
2. Triangulate the faces of Pᵢ into simplices.
3. Combine each boundary simplex with an interior point in Pᵢ. The set of
all such combinations form a triangulation of Pᵢ.
4. Calculate the volume of each simplex in the resulting triangulation. The
sum of these volumes is the volume of the intersection.
"""
function simplexintersection(S1::Array{Float64, 2}, S2::Array{Float64, 2}; tol::Float64 = 1/10^10)

  # Dimension
  const n = size(S1, 1)

  # Centroid and radii
  const c1 = Circumsphere(S1)[2:n+1]
  const c2 = Circumsphere(S2)[2:n+1]
  const r1 = Circumsphere(S1)[1]
  const r2 = Circumsphere(S2)[1]

  # Orientation of simplices
  const orientation_S1 = det([ones(1, n + 1); S1])
  const orientation_S2 = det([ones(1, n + 1); S2])

  if abs(orientation_S1) < tol || abs(orientation_S2) < tol
    return 0
  end

  # Set volume to zero initially. Change only if there is intersection
  IntVol = 0.0

# -------------------------------------
# Simplices intersect in some way
# -------------------------------------

  # If the (distance between centroids)^2-(sum of radii)^2 < 0,
  # then the simplices intersect in some way.
  dist_difference::Float64 = ((c1 - c2).' * (c1 - c2) - (r1 + r2)^2)[1]

  if dist_difference < 0
    # Find the number of points of each simplex contained within the
    # circumsphere of the other simplex

    vertices1InCircum2 = SomeVertexInCircumsphere(S1, r2, c2)
    vertices2InCircum1 = SomeVertexInCircumsphere(S2, r1, c1)

    # At least one circumsphere contains vertices of the other simplex
    #println("")

    if vertices1InCircum2 + vertices2InCircum1 >= 1
      #print("Finding barycentric coordinates\t")
      βs1in2, βs2in1, ordered_vertices1, ordered_vertices2, numof1in2, numof2in1 =
        BarycentricCoordinates(S1,S2,orientation_S1,orientation_S2,tol)
      # Trivial intersections
      TriviallyContained = heaviside0([numof1in2 numof2in1] - (n+1))
      IsSomeContained = sum(TriviallyContained, 2)[1]

      if IsSomeContained == 2.0 # The simplices coincide
        IntVol = abs(orientation_S1)
      elseif IsSomeContained == 1.0 # One simplex is contained in the other
        if TriviallyContained[1] == 1.0 # Simplex1 is contained in Simplex2
          IntVol = abs(orientation_S1)
        else # Simplex2 is contained in Simplex1
          IntVol = abs(orientation_S2)
        end
      else # No simplex contains the other

        #print("SharedVertices\t\t\t")

        Ncomm, ordered_vertices1, ordered_vertices2 = SharedVertices(βs1in2,ordered_vertices1,ordered_vertices2,numof1in2,numof2in1)
        # Is there any shared face?
        if Ncomm == n
          #print("The simplices share a face\t")
          IntVol = SharedFaceVolume(S2, βs1in2, ordered_vertices1, ordered_vertices2)
        else # The simplices do not share a face.
          #println("\nThe simplices do not share a face")
          #print("Intersection of boundaries\t")
          #@time IntVert, ConvexExpIntVert  = IntersectionOfBoundaries(S1,S2,βs1in2,βs2in1, ordered_vertices1, ordered_vertices2, numof1in2, numof2in1, Ncomm, tolerance)

          IntVert, ConvexExpIntVert = IntersectionOfBoundaries_NoStorage(S1,S2,βs1in2,βs2in1, ordered_vertices1, ordered_vertices2, numof1in2, numof2in1, Ncomm, tol)
          if !isempty(IntVert)
            #print("PolytopeGeneratingVertices\t")

            IntVert,ConvexExpIntVert = PolytopeGeneratingVertices(S1,S2,IntVert,ConvexExpIntVert,βs1in2,βs2in1,ordered_vertices1,ordered_vertices2,numof1in2,numof2in1,Ncomm);
            #print("Volume computation\t\t")
            IntVol = VolumeComputation(IntVert, ConvexExpIntVert)
            #println("")
          end
        end
      end
    else
      #println("No circumsphere of either simplex contains vertices of the other simplex")
    end
  end

  return IntVol
end


function intersectingvertices(S1::Array{Float64, 2}, S2::Array{Float64, 2}; tolerance::Float64 = 1/10^10, what = "volume")

  # Dimension
  const n = size(S1, 1)

  # Centroid and radii
  const c1 = Circumsphere(S1)[2:n+1]
  const c2 = Circumsphere(S2)[2:n+1]
  const r1 = Circumsphere(S1)[1]
  const r2 = Circumsphere(S2)[1]

  # Orientation of simplices
  const orientation_S1 = det([ones(1, n + 1); S1])
  const orientation_S2 = det([ones(1, n + 1); S2])

  if abs(orientation_S1) < tolerance || abs(orientation_S2) < tolerance
    return Float64[]
  end

  # If the (distance between centroids)^2-(sum of radii)^2 < 0,
  # then the simplices intersect in some way.
  dist_difference::Float64 = ((c1 - c2).' * (c1 - c2) - (r1 + r2)^2)[1]

  if dist_difference < 0
    # Find the number of points of each simplex contained within the
    # circumsphere of the other simplex

    vertices1InCircum2 = SomeVertexInCircumsphere(S1, r2, c2)
    vertices2InCircum1 = SomeVertexInCircumsphere(S2, r1, c1)

    # At least one circumsphere contains vertices of the other simplex
    #println("")

    if vertices1InCircum2 + vertices2InCircum1 >= 1
      #print("Finding barycentric coordinates\t")
      βs1in2, βs2in1, ordered_vertices1, ordered_vertices2, numof1in2, numof2in1 = BarycentricCoordinates(S1,S2,orientation_S1,orientation_S2,tolerance)

      # Trivial intersections
      TriviallyContained = heaviside0([numof1in2 numof2in1] - (n+1))
      IsSomeContained = sum(TriviallyContained, 2)[1]

      if IsSomeContained == 2.0 # The simplices coincide
        println("# The simplices coincide")
        return S1.'
      elseif IsSomeContained == 1.0 # One simplex is contained in the other
        if TriviallyContained[1] == 1.0 # Simplex1 is contained in Simplex2
          println("Simplex1 is contained in Simplex2")
          return S1.'
        else # Simplex2 is contained in Simplex1
          println("Simplex2 is contained in Simplex1")
          return S2.'
        end
      else # No simplex contains the other
        Ncomm, ordered_vertices1, ordered_vertices2 = SharedVertices(βs1in2,ordered_vertices1,ordered_vertices2,numof1in2,numof2in1)

        if Ncomm == n # The simplices share a face
          println("The simplices share a face")
          return SharedFaceVertices(S2, βs1in2, ordered_vertices1, ordered_vertices2)
        else  # The simplices do not share a face.
          IntVert, ConvexExpIntVert = IntersectionOfBoundaries_NoStorage(S1,S2,βs1in2,βs2in1, ordered_vertices1, ordered_vertices2, numof1in2, numof2in1, Ncomm, tolerance)
          if !isempty(IntVert)
            IntVert,ConvexExpIntVert = PolytopeGeneratingVertices(S1,S2,IntVert,ConvexExpIntVert,βs1in2,βs2in1,ordered_vertices1,ordered_vertices2,numof1in2,numof2in1,Ncomm)
            return IntVert
          else
            return Float64[]
          end
        end
      end
    else # No circumsphere of either simplex contains vertices of the other simplex
      return Float64[]
    end
  else
    return Float64[]
  end
  return Float64[]
end

function intersection(S1::Array{Float64, 2}, S2::Array{Float64, 2},
                      c1::Vector{Float64}, c2::Vector{Float64},
                      r1::Float64, r2::Float64,
                      orientation_S1::Float64, orientation_S2::Float64)
  tolerance = 1.0/10^10

  if abs(orientation_S1) < tolerance || abs(orientation_S2) < tolerance
    return 0
  end

  # Set volume to zero initially. Change only if there is intersection
  IntVol = 0.0

# -------------------------------------
# Simplices intersect in some way
# -------------------------------------

  # If the (distance between centroids)^2-(sum of radii)^2 < 0,
  # then the simplices intersect in some way.
  dist_difference::Float64 = ((c1 - c2).' * (c1 - c2) - (r1 + r2)^2)[1]

  if dist_difference < 0
    # Find the number of points of each simplex contained within the
    # circumsphere of the other simplex

    vertices1InCircum2 = SomeVertexInCircumsphere(S1, r2, c2)
    vertices2InCircum1 = SomeVertexInCircumsphere(S2, r1, c1)

    # At least one circumsphere contains vertices of the other simplex

    if vertices1InCircum2 + vertices2InCircum1 >= 1
      βs1in2, βs2in1, ordered_vertices1, ordered_vertices2, numof1in2, numof2in1 = BarycentricCoordinates(S1,S2,orientation_S1,orientation_S2,tolerance)

      # Trivial intersections
      TriviallyContained = heaviside0([numof1in2 numof2in1] - (n+1))
      IsSomeContained = sum(TriviallyContained, 2)[1]

      if IsSomeContained == 2.0 # The simplices coincide
        IntVol = abs(orientation_S1)
      elseif IsSomeContained == 1.0 # One simplex is contained in the other
        if TriviallyContained[1] == 1.0 # Simplex1 is contained in Simplex2
          IntVol = abs(orientation_S1)
        else # Simplex2 is contained in Simplex1
          IntVol = abs(orientation_S2)
        end
      else # No simplex contains the other

        @time Ncomm, ordered_vertices1, ordered_vertices2 = SharedVertices(βs1in2,ordered_vertices1,ordered_vertices2,numof1in2,numof2in1)


        # Is there any shared face?
        if Ncomm == n
          #print("The simplices share a face\t")
          @time IntVol = SharedFaceVolume(S2, βs1in2, ordered_vertices1, ordered_vertices2)
        else # The simplices do not share a face.
          #println("\nThe simplices do not share a face")
          #print("Intersection of boundaries\t")
          @time IntVert, ConvexExpIntVert  = IntersectionOfBoundaries(S1,S2,βs1in2,βs2in1, ordered_vertices1, ordered_vertices2, numof1in2, numof2in1, Ncomm, tolerance)
          if !isempty(IntVert)
            print("PolytopeGeneratingVertices\t")

            @time IntVert,ConvexExpIntVert = PolytopeGeneratingVertices(S1,S2,IntVert,ConvexExpIntVert,βs1in2,βs2in1,ordered_vertices1,ordered_vertices2,numof1in2,numof2in1,Ncomm);
            print("Volume computation\t\t")
            IntVol = VolumeComputation(IntVert, ConvexExpIntVert)
            println("")
          end
        end
      end
    else
      #println("No circumsphere of either simplex contains vertices of the other simplex")
    end
  end

  return IntVol
end


function simplexintersection(S1::Array{Float64, 2},
                            S2::Array{Float64, 2},
                            c1::Vector{Float64},
                            c2::Vector{Float64},
                            r1::Float64, r2::Float64,
                            orientation_S1::Float64,
                            orientation_S2::Float64)
  tolerance::Float64 = 1/10^10
  # Dimension
  const n = size(S1, 1)

  if abs(orientation_S1) < tolerance || abs(orientation_S2) < tolerance
    return 0
  end

  # Set volume to zero initially. Change only if there is intersection
  IntVol = 0.0

# -------------------------------------
# Simplices intersect in some way
# -------------------------------------

  # If the (distance between centroids)^2-(sum of radii)^2 < 0,
  # then the simplices intersect in some way.
  dist_difference::Float64 = ((c1 - c2).' * (c1 - c2) - (r1 + r2)^2)[1]

  if dist_difference < 0
    # Find the number of points of each simplex contained within the
    # circumsphere of the other simplex

    vertices1InCircum2 = SomeVertexInCircumsphere(S1, r2, c2)
    vertices2InCircum1 = SomeVertexInCircumsphere(S2, r1, c1)

    # At least one circumsphere contains vertices of the other simplex
    #println("")

    if vertices1InCircum2 + vertices2InCircum1 >= 1
      #print("Finding barycentric coordinates\t")
      βs1in2, βs2in1, ordered_vertices1, ordered_vertices2, numof1in2, numof2in1 = BarycentricCoordinates(S1,S2,orientation_S1,orientation_S2,tolerance)

      # Trivial intersections
      TriviallyContained = heaviside0([numof1in2 numof2in1] - (n+1))
      IsSomeContained = sum(TriviallyContained, 2)[1]

      if IsSomeContained == 2.0 # The simplices coincide
        IntVol = abs(orientation_S1)
      elseif IsSomeContained == 1.0 # One simplex is contained in the other
        if TriviallyContained[1] == 1.0 # Simplex1 is contained in Simplex2
          IntVol = abs(orientation_S1)
        else # Simplex2 is contained in Simplex1
          IntVol = abs(orientation_S2)
        end
      else # No simplex contains the other

        #print("SharedVertices\t\t\t")

        Ncomm, ordered_vertices1, ordered_vertices2 = SharedVertices(βs1in2,ordered_vertices1,ordered_vertices2,numof1in2,numof2in1)
        # Is there any shared face?
        if Ncomm == n
          #print("The simplices share a face\t")
          IntVol = SharedFaceVolume(S2, βs1in2, ordered_vertices1, ordered_vertices2)
        else # The simplices do not share a face.
          #println("\nThe simplices do not share a face")
          #print("Intersection of boundaries\t")
          #@time IntVert, ConvexExpIntVert  = IntersectionOfBoundaries(S1,S2,βs1in2,βs2in1, ordered_vertices1, ordered_vertices2, numof1in2, numof2in1, Ncomm, tolerance)

          IntVert, ConvexExpIntVert = IntersectionOfBoundaries_NoStorage(S1,S2,βs1in2,βs2in1, ordered_vertices1, ordered_vertices2, numof1in2, numof2in1, Ncomm, tolerance)
          if !isempty(IntVert)
            #print("PolytopeGeneratingVertices\t")

            IntVert,ConvexExpIntVert = PolytopeGeneratingVertices(S1,S2,IntVert,ConvexExpIntVert,βs1in2,βs2in1,ordered_vertices1,ordered_vertices2,numof1in2,numof2in1,Ncomm);
            #print("Volume computation\t\t")
            IntVol = VolumeComputation(IntVert, ConvexExpIntVert)
            #println("")
          end
        end
      end
    else
      #println("No circumsphere of either simplex contains vertices of the other simplex")
    end
  end

  return IntVol
end

function intersection2(S1::Array{Float64, 2}, S2::Array{Float64, 2},
                      c1::Vector{Float64}, c2::Vector{Float64},
                      r1::Float64, r2::Float64,
                      orientation_S1::Float64, orientation_S2::Float64)
  tolerance = 1.0/10^10

  if abs(orientation_S1) < tolerance || abs(orientation_S2) < tolerance
    return 0
  end

  # Set volume to zero initially. Change only if there is intersection
  IntVol = 0.0

# -------------------------------------
# Simplices intersect in some way
# -------------------------------------

  # If the (distance between centroids)^2-(sum of radii)^2 < 0,
  # then the simplices intersect in some way.
  dist_difference::Float64 = ((c1 - c2).' * (c1 - c2) - (r1 + r2)^2)[1]

  if dist_difference < 0
    # Find the number of points of each simplex contained within the
    # circumsphere of the other simplex

    vertices1InCircum2 = SomeVertexInCircumsphere(S1, r2, c2)
    vertices2InCircum1 = SomeVertexInCircumsphere(S2, r1, c1)

    # At least one circumsphere contains vertices of the other simplex

    if vertices1InCircum2 + vertices2InCircum1 >= 1
      βs1in2, βs2in1, ordered_vertices1, ordered_vertices2, numof1in2, numof2in1 = BarycentricCoordinates(S1,S2,orientation_S1,orientation_S2,tolerance)

      # Trivial intersections
      TriviallyContained = heaviside0([numof1in2 numof2in1] - (n+1))
      IsSomeContained = sum(TriviallyContained, 2)[1]

      if IsSomeContained == 2.0 # The simplices coincide
        IntVol = abs(orientation_S1)
      elseif IsSomeContained == 1.0 # One simplex is contained in the other
        if TriviallyContained[1] == 1.0 # Simplex1 is contained in Simplex2
          IntVol = abs(orientation_S1)
        else # Simplex2 is contained in Simplex1
          IntVol = abs(orientation_S2)
        end
      else # No simplex contains the other

        @time Ncomm, ordered_vertices1, ordered_vertices2 = SharedVertices(βs1in2,ordered_vertices1,ordered_vertices2,numof1in2,numof2in1)


        # Is there any shared face?
        if Ncomm == n
          #print("The simplices share a face\t")
          @time IntVol = SharedFaceVolume(S2, βs1in2, ordered_vertices1, ordered_vertices2)
        else # The simplices do not share a face.
          #println("\nThe simplices do not share a face")
          #print("Intersection of boundaries\t")
          @time IntVert, ConvexExpIntVert  = IntersectionOfBoundaries(S1,S2,βs1in2,βs2in1, ordered_vertices1, ordered_vertices2, numof1in2, numof2in1, Ncomm, tolerance)
          if !isempty(IntVert)
            print("PolytopeGeneratingVertices\t")

            @time IntVert,ConvexExpIntVert = PolytopeGeneratingVertices(S1,S2,IntVert,ConvexExpIntVert,βs1in2,βs2in1,ordered_vertices1,ordered_vertices2,numof1in2,numof2in1,Ncomm);
            print("Volume computation\t\t")
            IntVol = VolumeComputation(IntVert, ConvexExpIntVert)
            println("")
          end
        end
      end
    else
      #println("No circumsphere of either simplex contains vertices of the other simplex")
    end
  end

  return IntVol
end
