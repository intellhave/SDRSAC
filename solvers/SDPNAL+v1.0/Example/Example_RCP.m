%%*********************************************************************
%% Solve SDP relaxation of a clustering problem using the interface
%%
%% max { <W,X> | X*e = e, <I,X> = K, X psd, X >=0}
%% SDPNAL+: 
%% Copyright (c) 2017 by
%% Yancheng Yuan, Kim-Chuan Toh, Defeng Sun and Xinyuan Zhao
%%*********************************************************************

   clear all;
   rawdata = importdata('abalone-small.dat');
   W0 = rawdata.data;
   W = W0*W0';
   K = 5; 
   n = length(W); 
   e = ones(n,1); 
model = ccp_model('RCP');
    X = var_sdp(n,n);
    model.add_variable(X);
    model.maximize(inprod(W,X));
    model.add_affine_constraint(trace(X) == K);    
    %% add the constraints X*e = e
    options=1;
    if (options==1)
       for k=1:n
          ek = zeros(n,1); ek(k)=1;
          Ak = e*ek';
          model.add_affine_constraint(inprod(Ak,X) == 1);
       end
    else %% another way to specify the constraints
       At = sparse(n*n,n); 
       for k=1:n
           ek = sparse(n,1); ek(k)=1; Ak = e*ek';
           At(:,k) = Ak(:); 
       end
       model.add_affine_constraint(At'*map_vec(X)==e); 
    end
    model.add_affine_constraint(X >= 0);
model.solve;
Xsol = get_value(X);
%%*********************************************************************