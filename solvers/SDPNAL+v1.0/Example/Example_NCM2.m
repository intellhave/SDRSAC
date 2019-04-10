%%********************************************************************
%% This is a test example on using the interface
%% to solve the following Nearest Correlation Matrix problem 
%% min_X { norm(X-G,1) | diag(X) = ones(n,1), X psd }
%%
%% 
%% SDPNAL+: 
%% Copyright (c) 2017 by
%% Yancheng Yuan, Kim-Chuan Toh, Defeng Sun and Xinyuan Zhao
%%********************************************************************

clear all
rng('default')
n = 100;
G = randn(n,n);
G = 0.5*(G' + G);
H = rand(n,n);
blk{1,1} = 's'; blk{1,2} = n; 
n2 = n*(n+1)/2; 
model = ccp_model('Example_NCM2');
    X  = var_sdp(n,n);
    Xp = var_symm(n,n);
    Xn = var_symm(n,n);
    model.add_variable(X,Xp,Xn);
    model.minimize(inprod(ones(n),Xp)+inprod(ones(n),Xn));
    model.add_affine_constraint(map_diag(X) == ones(n,1));
    %%model.add_affine_constraint(X-Xp+Xn==G);
    model.add_affine_constraint(map_svec(X)-map_svec(Xp)+map_svec(Xn)==svec(blk,G));
    model.add_affine_constraint(Xp >= 0); 
    model.add_affine_constraint(Xn >= 0); 
model.solve;
%%********************************************************************