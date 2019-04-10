%%********************************************************************
%% This script is for QAP
%% SDPNAL+: 
%% Copyright (c) 2017 by
%% Yancheng Yuan, Kim-Chuan Toh, Defeng Sun and Xinyuan Zhao
%%********************************************************************

clear all;
problem_name = 'chr20b';
[A, B] = qapread(strcat(problem_name, '.dat'));
%% Construct C
Ascale = max(1, norm(A, 'fro'));
Bscale = max(1, norm(B, 'fro'));
A = A/Ascale; B = B/Bscale;
C = kron(B, A); C = 0.5*(C + C');
nn = length(C);
n = length(A);
%%===================================================================
model = ccp_model(['QAP-',problem_name]);
Y = var_sdp(nn, nn);
model.add_variable(Y);
model.minimize(inprod(C, Y));
model.add_affine_constraint(Y >= 0);
II = speye(n,n); 
EE = ones(n,n); 
for i = 1:n-1
    for j = i:n
        Eij = sparse(i,j,1,n,n); 
        if i == j
           const = 1;
        else
           const = 0;
        end
        model.add_affine_constraint(inprod(kron(II,Eij), Y) == const);
        model.add_affine_constraint(inprod(kron(Eij,II), Y) == const);
        model.add_affine_constraint(inprod(kron(Eij,EE), Y) == 1);        
    end
end
model.add_affine_constraint(inprod(kron(II,sparse(n,n,1,n,n)), Y) == 1);
model.setparameter('maxiter',6000,'printlevel',2,'stopoption',0);
model.solve;
%%********************************************************************        