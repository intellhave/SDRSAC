%%***********************************************************************
%% NCM: nearest correlation matrix
%% min { sum(W.*abs(X-G)) : diag(X) = e, X psd}
%%
%% X-G = U-V, with U,V >= 0.
%%
%% W = a nonnegative weight matrix. 
%%***********************************************************************
%% SDPNAL+ 
%% Copyright (c) 2014 by
%% Liuqin Yang, Defeng Sun, and Kim-Chuan Toh
%%***********************************************************************

   function [blk,At,C,b] = NCM_Hnorm(G,W,options)

   n = length(G);
   if (nargin < 3); options = 1; end
   if (nargin < 2); W = ones(n); end

   n2 = n*(n+1)/2; 
   blk{1,1} = 's'; blk{1,2} = n; 

   ww = svec(blk(1,:),W); 
   gg = svec(blk(1,:),G); 

   C{1,1} = sparse(n,n);    
   Acell = cell(1,n); 
   for k = 1:n; Acell{k} = spconvert([k,k,1;n,n,0]); end
   Atmp = svec(blk(1,:),Acell,1); 

   if (options==1)
      II = speye(n2,n2); 
      blk{2,1} = 'l'; blk{2,2} = 2*n2; 
      At{1,1} = [Atmp{1}, II];       
      At{2,1} = [sparse(n,2*n2); -II, II]'; 
      C{2,1} = [ww; ww];       
      b = [ones(n,1); gg];       
   else
      %% removed the diagonal constraints
      %% in X-U+V = G
      tmp = svec(blk(1,:),speye(n,n)); 
      idx = find(tmp==0); 
      II = speye(n2,n2); 
      II = II(idx,:);
      blk{2,1} = 'l'; blk{2,2} = 2*n2-2*n; 
      At{1,1} = [Atmp{1}, II'];       
      At{2,1} = [sparse(n,2*n2-2*n); -speye(n2-n) speye(n2-n)]';       
      C{2,1} = [ww(idx); ww(idx)]; 
      b = [ones(n,1); gg(idx)]; 
   end
%%*****************************************************
