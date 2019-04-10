function inls = countCorrespondences(M, Btree, epsilon)

    [idx, dist] = knnsearch(Btree, M');
    in_idx = find(dist<=epsilon);
    
    inls = length(idx(in_idx));

end