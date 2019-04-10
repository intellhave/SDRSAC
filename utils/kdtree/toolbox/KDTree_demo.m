clc, clear, close all;
 
% Input data
data = rand(1000,2);
 
% Build tree (Memory automatically freed)
kd = KDTree(data);
 
% Nearest neighbors
query = rand(10,2);
nn_idxs = kd.nn(query);
% disp(data(nn_idxs,:)); 
 
% K nearest neighbors
K = 10;
query = rand(1,2);
knn_idxs = kd.knn(query,K);
% disp(data(knn_idxs,:)); 

% Ball query
radii = .5;
query = rand(1,2);
ball_idxs = kd.ball(query, radii);
% disp(data(ball_idxs,:));

% Support for save
save matlab.mat kd;

% Now simply reload without re-building!!
clear all;
load matlab.mat kd;
ball_idxs = kd.ball(rand(1,2), .5);

delete matlab.mat