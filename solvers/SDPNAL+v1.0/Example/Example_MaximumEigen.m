%%*********************************************************************
%% This is an example to test the Maximum Eigenvalue problem.
%% min_y { max(eig(C + Ly)) }
%% s.t. norm(y,inf) <= 1
%% SDPNAL+: 
%% Copyright (c) 2017 by
%% Yancheng Yuan, Kim-Chuan Toh, Defeng Sun and Xinyuan Zhao
%%*********************************************************************

clear all;
rng('default')
n = 25;
p = 10*n^2;

C = randn(n,n);
C = 0.5*(C + C');
LT = cell(p,1);
for i = 1:p
    LT{i} = randn(n,n);
    LT{i} = 0.5*(LT{i}+LT{i}');
end

model = ccp_model('Maximum_Eigenvalue');
    X = var_sdp(n,n);
    u = var_nn(p);
    v = var_nn(p);
    model.add_variable(X, u, v);
    model.minimize(inprod(-C, X) + inprod(u,ones(p,1)) + inprod(v,ones(p,1)));
    model.add_affine_constraint(LT*X - u + v == 0);
    model.add_affine_constraint(trace(X)==1);
    model.setparameter('AATsolve.method','iterative', 'maxiter', 2000);
    model.setparameter('printlevel',2);
model.solve
XX = get_value(X); 
uu = get_value(u);
vv = get_value(v); 
yy = model.info.dual_opt.y(1:end-1);
%*********************************************************************
if (false)
   blk=model.info.prob.blk; 
   At=model.info.prob.At; C=model.info.prob.C; b=model.info.prob.b;
   L=model.info.prob.L; U=model.info.prob.U;
end
%*********************************************************************