% Given two point set M and B, find the set of correspondences using SDP

function [R,t, corrM, corrB] = sdpReg(M, B, config)
    
    % Generate Weight Matrix
    corrM = []; corrB = []; R=[]; t = [];    
    
    [W, acceptedPairs] = generateWeightMatrix(M, B, config);    
           
    if (acceptedPairs > config.nPairThresh)
        
        k = config.k;
        
        [X] = solveSDP(W, k);    
        
        [corrM, corrB] = getCorrespondences(X);
        sM = M(:,corrM);
        sB = B(:, corrB);    
        [R, t] = Kabsch(sM, sB);
        
    end
    
    

end

