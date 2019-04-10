% Main SDRSAC Algorithm
% Input: 
%       M, B: Input point clouds, each is a 3xN matrix  
%       config: Contains all parameters required to run the algoirthm.
%               See readConfig.m for more info    
% Output: a structure out containing output variables:
%         out.bestR: Best rotation matrix
%         out.bestT: Best translation vector


function [out] = SDRSAC(M, B, config)    
    add_dependencies;    
    ps = 0.99;
            
    iter=0;
    T_max = 1e10;
    
    % Using KDTree for quick computing of consensus size
    B_Tree = KDTreeSearcher(B');
    
    % Prepare sampling
    n = config.pointPerSample;  % N_sample 
    maxInls = 0;
    bestR = [];
    bestT = [];
    
    % Start the sampling iterations
    stop = false;
    while (iter < config.maxIter && ~stop)        
        
        idxM = randsample(size(M,2), n);
        m = M(:, idxM);       
        B_to_sample = B;
        scount = 0;
        while (size(B_to_sample,2) > n*2 && scount < 2)
            scount = scount + 1;
            fprintf('Current B size: %d\n', size(B_to_sample,2));
            
            idxB = randsample(size(B_to_sample,2), n);                    
            b = B_to_sample(:, idxB);                        
            B_to_sample(:,idxB) = [];
            
            % Solve SDP
            [Rs,ts, ~, corrB] = sdpReg(m, b, config);

            if ~isempty(corrB)
                TM = Rs*M + repmat(ts, 1, size(M,2));
               
                
                % Conducting ICP
                [Ricp, Ticp] = icp(B, TM, 'Matching', 'kDtree', 'WorstRejection', 0.1, 'iter', 100);               
                TMICP = Ricp*Rs*M + repmat(Ricp*ts + Ticp, 1, size(M, 2));                                                               
                inls_icp = countCorrespondences(TMICP, B_Tree, config.epsilon);   
                
                if inls_icp > maxInls
                    maxInls = inls_icp;
                    bestR = Ricp*Rs;
                    bestT = Ricp*ts + Ticp;  
                                       
                    fprintf('Best-so-far consensus size: %d\n', maxInls);%
                    fprintf('--------------');
                    
                    % For debugging purpose:
                     %close all; plotPointClouds(B, TMICP, 'b.','r.');
                    
                    % Compute stopping criterion
                    pI = maxInls./size(M,2);
                    T_max = log(1-ps)./log(1-pI^config.k);
                end                
            end 
            
            iter = iter+1;
            
            % Stopping criterion is satisfied (minimum of 5 iterations)
            if iter >= T_max && iter >=5 
                 stop = true;
            end
        end
        
    end
    
    out.inls = maxInls;    
    out.iter = iter;    
    out.R = bestR;
    out.T = bestT;    
    plotPointClouds(B, bestR*M + bestT, 'b.','r.')
end 
    