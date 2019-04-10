%%***********************************************************************
%% biq: generate SDP data for binary quadratic programming
%%
%% X = [Y,y; y', 1];
%%***********************************************************************
%% SDPNAL+ 
%% Copyright (c) 2014 by
%% Liuqin Yang, Defeng Sun, and Kim-Chuan Toh
%%***********************************************************************

   function [blk,At,C,b] = biq(Q,options)

   if (nargin == 1); options = 0; end;
   if isempty(Q)
      blk = []; At = []; C = []; b = [];
      return;
   end
   if norm(Q-Q','fro') 
      Q = 0.5*(Q + Q'); 
   end
   
   n = size(Q,1);   
   blk{1,1} = 's'; blk{1,2} = n+1;   
%%
%% 
%%
   C{1} = [Q,sparse(n,1); sparse(1,n+1)];  
%%
%% diag(Y) = y; 
%% 
   n1 = n+1;   
   Acell = cell(1,n1); 
   for k = 1:n
      Acell{k} = spconvert([k,k,1;k,n1,-0.5; n1,k,-0.5;n1,n1,0]); 
   end
   Acell{n1} = spconvert([n1,n1,1]); 
   At = svec(blk,Acell,1); 
   b = [zeros(n,1); 1]; 
%%
%% add the constraint <E,Y> >= (e'*y)^2
%%
   if (options) %% Not useful. 
      e = ones(n,1); 
      tmp{1,1} = [sparse(n,n), 0.5*e; 0.5*e', 0];  
      tmp{1,2} = [ones(n,n), zeros(n,1); zeros(1,n), 0];  
      tmp{1,3} = sparse(n+1,n+1); 
      Atmp = svec(blk(1,:),tmp,1); 
      At{1} = [At{1}, Atmp{1}]; 
      b = [b; 0;0;1];
      C{2,1} = zeros(3,1); 
      blk{2,1} = 'q'; blk{2,2} = 3; 
      At{2,1} = [sparse(n1,3); 0,0,-1; -1,-1,0; 1,-1,0]'; 
   end
%%******************************************************
