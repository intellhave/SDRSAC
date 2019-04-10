function syntheticExperiment(method)

if (nargin<1)
    method = 'SDRSAC';
end

method = 'SDRSAC';
%method = '4PCS';
% Synthetic experiment

config = readConfig();
%config = readConfig('redwood');

load(config.matPath);

%outRANSAC = pointCloudReg(D, M, config, 'RANSAC');
%outICP = pointCloudReg(M, D, config, 'ICP');

% if (strcmp(method,'SDRSAC'))
% 
%     D = normalize(D,2, 'center');
%     M = normalize(M,2, 'center');
%     % 
%     D = normalize(D,2, 'range');
%     M = normalize(M,2, 'range');
% end

%out = pointCloudReg(D, M, config, method);    
%plotPointClouds(D, M, 'r.','b.');
log_file_path = './logs/logs.txt';
for nrun = 1:10
    for m = 3
        config.k = m;
        for nsample = 16
            config.pointPerSample = nsample;
            out = pointCloudReg(D, M, config, method);
            disp(out);        
            if strcmp(method,'SDRSAC')
                log_results(log_file_path, out, config);
            end
        end
    end
end

disp(out);

end




%outSDR = pointCloudReg(M, D, config, '4PCS');
%outGlobGM = pointCloudReg(M, D, config, 'GlobGM');
%[out] = SDRReg(D, M, config);