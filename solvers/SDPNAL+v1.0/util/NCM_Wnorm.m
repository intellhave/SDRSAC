%%**********************************************************************
%% NCM_Wnorm: solve the W-norm weighted nearest correlation 
%%            matrix problem. 
%%
%% Given correlation data matrix B, find X that solves: 
%%
%% min { sum abs(W*(X-B)*W) : diag(X)=1, X psd } 
%% ==
%% min { W(Y1+Y2)W : diag(X)=1, X-Y1+Y2=B, X psd } 
%% 
%% Input: 
%% B = data matrix, or 
%%     the name of the MATLAB *.mat file containing the data matrix B.
%% W = The weight matrix (optional). The default is W = I. 
%%***********************************************************************
%% SDPNAL+ 
%% Copyright (c) 2014 by
%% Liuqin Yang, Defeng Sun, and Kim-Chuan Toh
%%***********************************************************************

     function [blk,At,C,b] = NCM_Wnorm(B,W)

     if isstr(B); 
        eval(['load ',B]); 
     end
     n = length(B); 
     if (nargin == 1); W = speye(n,n); end
     if (norm(B-B','fro') > 1e-13)
        fprintf('warning: B is not symmetric; proceed to symmetrize it');
        B = 0.5*(B+B'); 
     end
     if (norm(W-W','fro') > 1e-13)
        fprintf('warning: W is not symmetric; proceed to symmetrize it');
        W = 0.5*(W+W'); 
     end
     [Wchol,indef] = chol(W+1e-8*speye(n));
     if (indef); error('W is not positive semidefinite'); end 
     dW = full(W*ones(n,1));      
     %%----------------------------------------------
     %% generate At,C,b 
     %%----------------------------------------------
     n2 = n*(n+1)/2; 
     blk{1,1} = 's'; blk{1,2} = n;
     blk{2,1} = 'l'; blk{2,2} = 2*n2;    
    
     Acell = cell(1,n);
     for k = 1:n;
        Acell{k} = spconvert([k,k,1; n,n 0]);
     end      
     Atmp = svec(blk(1,:),Acell,1); 

     II = speye(n2); 
     At{1,1} = [Atmp{1}, II]; 
     At{2,1} = [sparse(n,2*n2); -II, II]'; 

     b = [ones(n,1); svec(blk(1,:),B)];    
     ee = svec(blk(1,:),dW*dW'); 
     C{1,1} = sparse(n,n); 
     C{2,1} = [ee;ee]; 
%%**********************************************************************
