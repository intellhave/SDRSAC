%%*********************************************************************
%% Solve NCM using interface
%% to solve the following NCM problem
%% min ||X - G||_1
%% s.t. diag(X)==1 ; 
%%      X \in S^+_n
%% SDPNAL+: 
%% Copyright (c) 2017 by
%% Yancheng Yuan, Kim-Chuan Toh, Defeng Sun and Xinyuan Zhao
%%*********************************************************************

clear all;
rng('default')
n = 100;
G = randn(n,n);
G = 0.5*(G+G');
model = ccp_model('NCM');
    X = var_sdp(n,n);
    model.add_variable(X);
    model.minimize(l1_norm(X-G));
    model.add_affine_constraint(map_diag(X) == ones(n,1));
    %%model.add_psd_constraint( X >= 1e-3*speye(n));
model.solve;
%%*********************************************************************