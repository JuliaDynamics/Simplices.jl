function speed_test_nontrivial(dim,N)
    for i = 1:N
        S1,S2 = nontrivially_intersecting_simplices(dim)
        simplexintersection(S1.', S2.')
    end
end


function speed_test_sharing(dim,N)
    for i = 1:N
        S1,S2 = simplices_sharing_vertices(dim)
        simplexintersection(S1.', S2.')
    end
end
