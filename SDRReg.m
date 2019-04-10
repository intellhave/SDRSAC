%SDRSAC 
% M and B are input point cloud (can be dense)
function [out] = SDRReg(M, B, config)    
    add_dependencies;    
    %rng(2912670);       
    iter=0;
    stop = false;
    out.time = 0;  out.iter = 0;  out.maxInls = 0;       
    
    M_subsampled = subSamplePointCloud(M, config.NSubsample, 'random');
    B_subsampled = subSamplePointCloud(B, config.NSubsample, 'random');    
    B_subsampled_Tree=KDTreeSearcher(B_subsampled');
    %Now, compute the surface normals    
    normalM = computeNormal(M);    normalB = computeNormal(B);        
    MTree = KDTreeSearcher(M'); BTree = KDTreeSearcher(B');
    cMIdx = knnsearch(MTree, M_subsampled'); 
    cBIdx = knnsearch(BTree, B_subsampled');
    normalM_subsampled = normalM(:, cMIdx);
    normalB_subsampled = normalB(:, cBIdx);

    n = config.pointPerSample;
    maxInls = 0;

    % For large-span sampling
    bb = boundingBox(M_subsampled);
    maxDist = abs(bb(:,2) - bb(:,1));
    config.pointDistance = config.pointDistance*min(maxDist);             
    
    % For parallel   
    PTM = cell( config.innerIter,1); PRT = cell( config.innerIter, 1); 
    PInls = zeros( config.innerIter,1);   PTime = zeros(config.innerIter, 1);

    % Start sampling
    while (iter < config.maxIter&& ~stop)        
        tic;        
        idxBPool = []; idxM = [];                
        for sp = 1:config.kSample       
            CS = [];            
            while (isempty(CS))
                disp('---------Subsample for congruent----------');
                [idxMC, r1, r2, d1, d2, angle] = sampleCoplanarPoints(M_subsampled, idxM, config);
                B_subsampled_for_congruent = B_subsampled(:, randsample(size(B_subsampled,2), min(config.SamplesForCongruent, size(B_subsampled,2))));                    
                CS =  findCongruentSet(B_subsampled_for_congruent, r1, r2, d1, d2, angle, config);
                
%                 close all; plotPointClouds(M_subsampled, M_subsampled(:, idxMC), 'b.', 'ro');
%                 hold on;
%                 plot_plane( M_subsampled(:, idxMC(1:3)));
%                 if (length(CS)>0)
%                     figure;
%                     plotPointClouds(B_subsampled, B_subsampled(:, CS(1,:)), 'b.', 'ro');
%                     hold on;
%                     plot_plane( B_subsampled(:, CS(1,1:3)));
%                 end
%                 
            end                                          
            idxM = [idxM; idxMC];
            idxBPool = [idxBPool; CS(:)];
        end
        % Now, a subset of M has been sampled and its congruent sets
        idxBPool = unique(idxBPool);
        M_congruent = M_subsampled(:, idxM); 
        normalM_congruent = normalM_subsampled(:, idxM);
        B_congruent = B_subsampled_for_congruent(:, idxBPool);
        normalB_congruent = normalB_subsampled(:, idxBPool);                       
        m = M_congruent;
        normalm = normalM_congruent;
        nIter = config.innerIter;
        
        % Generate samples for parfor
        gIdx =zeros(nIter, n);
        for it=1:nIter
            gIdx(it,:) = randsample(size(B_congruent,2), n)';
            %gIdx(it,:) = samplePoints(B_congruent, n, [], config);
        end        
        if (size(idxBPool)<=n+1) 
            nIter = 1; 
            disp('idxBPool'); disp(size(idxBPool));
        end
        out.time = out.time+toc;
        iter = iter + nIter;
        

        parfor it = 1:nIter     
            tic;
            idxB = gIdx(it,:);
            b = B_congruent(:, idxB);  
            normalb = normalB_congruent(:, idxB);
            [Rs,ts, ~, corrB] = sdpReg(m, b, normalm, normalb, config);            
            if (~isempty(corrB))
                TM = Rs * M_subsampled + repmat(ts, 1, size(M_subsampled, 2));                                                                             
                %[Ricp, Ticp] = icp(B_subsampled, TM, 'Matching', 'kDtree', 'iter', 20);
                %Ricp = eye(3); Ticp = zeros(3,1);
                [Ricp, Ticp] = icp(B_subsampled, TM, 'Matching', 'kDtree', 'WorstRejection', 0.3, 'iter', 20);
                %[Ricp, Ticp] = huber_icp(B, TM, 20, 5, 1);
                %[Ricp, Ticp] = huber_icp(B_subsampled, TM);
                %Ricp = Rs; Ticp = ts;                              
                TMICP = Ricp*Rs*M_subsampled + repmat(Ricp*ts + Ticp, 1, size(M_subsampled, 2));                
                inls_icp = countCorrespondences(TMICP, B_subsampled_Tree, config.epsilon_sampled);                 
                TMTree = KDTreeSearcher(TMICP');
                inls_icp2 = countCorrespondences(B_subsampled, TMTree,config.epsilon_sampled);
                inls_icp = min(inls_icp, inls_icp2);

                PRT{it} = [Ricp*Rs Ricp*ts + Ticp];                 
                %PRT{it} = [Rs ts];                 
                PTM{it} =  Ricp*Rs*M_subsampled + repmat(Ricp*ts + Ticp, 1, size(M_subsampled, 2));        
                PInls(it) = inls_icp;                                    
            end
            PTime(it) = toc;
            disp(PTime(it));
        end
        out.time = out.time + sum(PTime);
        disp(PInls); 
        disp(['======CURRENT MAX  = ' num2str(maxInls)]);
        [nInls,idx] = max(PInls);
        if (nInls > maxInls)                                  
            maxInls = nInls;                                              
            bestRT = PRT{idx}; 
            out.bestR = bestRT(:,1:3); out.bestT = bestRT(:, 4);
            TM = out.bestR*M + repmat(out.bestT, 1, size(M, 2));
            %[Ricp, Ticp] = icp(B, TM, 'Matching', 'kDtree', 'WorstRejection', 0.2, 'iter', 50);
            %[Ricp, Ticp] = huber_icp(B, TM);
%             [Ricp, Ticp] = icp(B, TM);
%             TMICP = Ricp * TM + Ticp;            
%             
            TMTree = KDTreeSearcher(TM');
            totalInls = countCorrespondences(TM, BTree, config.epsilon);    
            totalInls2 = countCorrespondences(B, TMTree, config.epsilon);    
            out.inls = maxInls;
            out.totalIns = totalInls;          
            close all; plotPointClouds(B_subsampled, PTM{idx}, 'b.','r.');  
            figure;
            plotPointClouds(B, TM, 'b.','r.');            
            disp('---MaxInls----'); disp(maxInls);     
            disp('---Total----'); disp(totalInls);      
            disp('---Total2----'); disp(totalInls2);      
            disp(out.bestR); disp(out.bestT);           
            out.cloudA =  TM; out.cloudB = B;            
        end
        T = computeStopping(0.99, maxInls./size(M,2));
        if (iter>=T && iter >= config.minIterRequired)
            stop = true; 
        end
    end
        
        
end

   