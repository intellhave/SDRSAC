%%*******************************************************************
%% randSDP.m : creates random feasible SDP problems of the form:
%%
%% min { <C,X> : A(X) = b, X psd} 
%% where A,C,b are random data.
%%
%% [blk,At,C,b] = randSDP(n,m);
%% n = matrix dimension
%% m = number of linear constraints, i.e. dimension of b. 
%%***********************************************************************
%% SDPNAL+ 
%% Copyright (c) 2014 by
%% Liuqin Yang, Defeng Sun, and Kim-Chuan Toh
%%***********************************************************************

   function  [blk,At,C,b] = randSDP(n,m)

   randn('state',0); 
   rand('state',0); 

   blk{1,1} = 's'; blk{1,2} = n; 
%%
%% generate feasible point X0,Z0 
%% 
   tmp = randn(n); tmp = tmp*tmp';
   X0 = 0.5*(tmp + tmp');
   tmp = randn(n); tmp = tmp*tmp';
   Z0 = 0.5*(tmp + tmp');
%%
%% set up the matrices Ak,C,b
%%
   b = zeros(m,1); 
   y0 = randn(m,1); 
   AA  = cell(1,m);
   Aty0 = zeros(n,n);   
   for k = 1:m;
      tmp = randn(n,n); Ak = 0.5*(tmp+tmp');
      AA{1,k} = Ak;
      b(k) = sum(sum(Ak.*X0));  
      Aty0 = Aty0 + y0(k)*Ak; 
   end 
   C{1} = Z0 + Aty0; 
   At = svec(blk,AA); 
%%******************************************************************
