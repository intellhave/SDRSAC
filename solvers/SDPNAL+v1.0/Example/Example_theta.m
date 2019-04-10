%%*********************************************************************
%% This is an example of using the interface
%% to solve the following theta problem of a graph 
%% with adjacency matrix G: 
%% max <E,X>
%% s.t. trace(X)==1 ; X positive semidefinite
%%      X(i,j) == 0 for (i,j) in Edge_set
%% SDPNAL+: 
%% Copyright (c) 2017 by
%% Yancheng Yuan, Kim-Chuan Toh, Defeng Sun and Xinyuan Zhao
%%*********************************************************************

 load theta6.mat
 [IE,JE] = find(triu(G,1));
 n = length(G); 
 %%
 model = ccp_model('Example_theta');
     X = var_sdp(n, n);
     model.add_variable(X);
     model.maximize(sum(X));
     model.add_affine_constraint(trace(X) == 1);
     options=2; 
     if (options==1)
        model.add_affine_constraint(X(IE,JE) == 0);
        model.add_affine_constraint(X >= 0);
     else %% another way to specify the constraints
        U = 1e20*(1-spones(G+G'));        
        model.add_affine_constraint(0 <= X <= U); 
     end
 model.solve;
%%*********************************************************************