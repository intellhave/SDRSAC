% Get corrspondences from a permutation matrix X
% X is a permuation matrix

function [corrM, corrB] = getCorrespondences(X)
   
    %n = size(X, 1);
    [y, idx] = max(X, [], 2);
    corrM = find(y>0);
    corrB = idx(corrM);
    %[~,corrB] =  max(X, [], 2);    
    
   
    
end