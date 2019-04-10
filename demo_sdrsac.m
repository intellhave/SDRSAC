% Demo Program for the paper: 
%   SDRSAC: Semidefinite-Based Randomized Approach for Robust Point Cloud
%   Registration without Correspondences - Huu Le, Thanh-Toan Do, Tuan
%   Hoang and Ngai-Man Cheung - CVPR 2019
%       
%   
dataset = 'synthetic';
method = 'SDRSAC';          % Without Correspondences


% Read configuration containing hyperparameters
config = readConfig(dataset); 

% Read Data
load(config.matPath);

% Run SDRSAC
out = pointCloudReg(D, M, config, method);

disp(out);
