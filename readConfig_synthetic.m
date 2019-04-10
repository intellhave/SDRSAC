function config = readConfig_synthetic()

    
   
    config.datasetName = 'bunny10';
    
    config.maxIter = 10000;
    config.innerIter = 10;
    
    config.k = 4;               % m in the paper
    config.pointPerSample = 16; % Nsample
    
    config.epsilon = 0.015;
    
    config.nPairThresh = 5;   % Sample is accepted if the number of pairs that has distance > dDiffThreshold is greater than this
    config.dDiffThresh = 0.0001; % Used for generating weight matrix
    config.pairDistThresh = 0.01;

         
    % For gsynthetic data   
    config.syntheticN = 5000;
    config.OutlierRate = 30;

    
    config.SamplesForCongruent = 500;
    config.extractPairEpsilon = 0.01;  % Grid size for extraing pairs
    
    config.angleDiffThreshold = 0.01; % For congruent testing
    
    %-------------DATASET-------------
    plyMaps = containers.Map(); plyMapsB = containers.Map();
    matMaps = containers.Map();
%     plyOutMaps = containers.Map(); 
    
    % Synthetic data with 10% outliers
    plyMaps('bunny10') = './data/bun045.ply';
    plyMapsB('bunny10') = './data/bun045.ply';
    matMaps('bunny10') = './data/bunny10.mat';
   % plyOutMaps('bunny10') = '../Dataset/bunny/data/bunny10_synthetic';
    %%% More datasets can be addded in a similar manner
    
    config.plyPath = plyMaps(config.datasetName);
    config.plyPathB = plyMapsB(config.datasetName);
    config.matPath = matMaps(config.datasetName);
end