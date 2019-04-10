% Generate weight matrix for the task of quadratic matching
function  [W, acceptedPairs] = generateWeightMatrix(M, B, config)
    
    N = size(M, 2);    
    W = -10000*ones(N*N, N*N);
    acceptedPairs = 0;
    
    % Loop through all possible pairs
    for p=1:N
        for q=1:N
            for s=1:N
                for t=1:N
                    if (p~=s && q~=t)
                        i = p + (q-1)*N;
                        j = s + (t-1)*N;
                        dDiff = abs( distance(M,p,s) - distance(B,q,t));                         
                        psDist = distance(M, p,s);
                        qtDist = distance(B, q,t);
                        
                        if dDiff <= config.dDiffThresh && psDist >= config.pairDistThresh && qtDist >= config.pairDistThresh                                     
                            W (i,j) = exp(-dDiff);  % can also use 1./dDiff; %
                            acceptedPairs = acceptedPairs + 1;                            
                        end                                                        
                    end                    
                end
            end
        end
    end
    

end