% Read point cloud from PLY file and return coordinates
function [X, Y, Z] = readPointCloud(plyFile, s)
    if (nargin < 2) s = 1; end;
    
    
    ptCloud = pcread(plyFile);
    N = ptCloud.Count;
    if (s<=1)
        idx = randsample(N, round(N*s))';
    else 
        idx = randsample(N, s)';
    end
    X = ptCloud.Location(idx,1)';
    Y = ptCloud.Location(idx,2)';
    Z = ptCloud.Location(idx,3)'; 
    
    X = double(X); Y = double(Y); Z = double(Z);
end