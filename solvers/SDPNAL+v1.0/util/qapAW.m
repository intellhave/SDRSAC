%%******************************************************
%% qapAW: generate SDP data for the problem QAP-AW1 in 
%%        the paper by Povh and Rendl.
%%
%% SDPNAL: 
%% Copyright (c) 2008 by
%% Xinyuan Zhao, Defeng Sun, and Kim-Chuan Toh 
%%******************************************************

   function [blk,Ft,CC,bb,Ascale,Bscale] = qapAW(A,B,options)
   
   n = size(A,1);
   if (size(B,1)~= n); 
      error('A and B must have the same size'); 
   end
   if (nargin == 2); options = 2; end
   
   n2 = n*n;
   blk{1,1} = 's'; blk{1,2} = n2;
%%
%% 
%%
   Ascale = max(1,norm(A,'fro')); 
   Bscale = max(1,norm(B,'fro'));
   A = A/Ascale; B = B/Bscale;
   CC{1} = kron(B,A); 
   if norm(A-A','fro') | norm(B-B','fro')
      CC{1} = 0.5*(CC{1} + CC{1}'); 
   end
%%
%% sum_i Y^{ii} = I_n
%% create n*(n+1)/2 constraints
%%   
   pp = [0:n:n*(n-1)]'; vv = (1/sqrt(2))*ones(n,1); 
   Fcell = cell(1,n*(n+1)/2);
   count = 0; 
   for j = 1:n
      for i = 1:j
         row = pp + i; col = pp + j; 
         if (i==j)
            Fcell{count+1} = spconvert([row,col,ones(n,1); n2,n2, 0]); 
         else
            Fcell{count+1} = spconvert([row,col,vv; n2,n2, 0]); 
         end
         count = count+1; 
      end
   end   
   Ft = svec(blk,Fcell,1); 
   bb = svec(blk,speye(n,n)); 
%%
%% <I_n,Y^{ij}> = delta_{ij} for j=1:n, i=1:j. 
%%   
   pp = [1:n]'; vv = ones(n,1)/2; 
   count = 0; 
   btmp = zeros(n*(n+1)/2,1); 
   for j = 1:n
      for i = 1:j
         row = pp + (i-1)*n; col = pp + (j-1)*n; 
         if (i==j) 
            Fcell{count+1} = spconvert([row,col,ones(n,1); n2,n2, 0]); 
            btmp(count+1) = 1; 
            count = count+1; 
         elseif (i < j)
            Fcell{count+1} = spconvert([row,col,vv; n2,n2, 0]); 
            btmp(count+1) = 0; 
            count = count+1; 
         end
      end
   end
   if (options == 1)
      %%
      %% add <E,Y> = n2
      %%
      idx = [1:n*(n+1)/2-1]; 
      alpha = n2; 
      ee = svec(blk,ones(n2,n2),1); 
      Ftmp = svec(blk,Fcell(idx),1); 
      Ft{1} = [Ft{1}, Ftmp{1}, ee/alpha];
      bb = [bb; btmp(idx); n2/alpha];  
   else
      %%
      %% add <ee^T,Y^{ij}> = 1 for j=1:n, i=1:j. 
      %% 
      Fcell2 = cell(1,n*(n+1)/2);
      alpha = sqrt(n);
      btmp2 = ones(n*(n+1)/2,1)/alpha;      
      rr = [1:n]'*ones(1,n); rr = rr(:);  
      cc = ones(n,1)*[1:n];  cc = cc(:);
      vv = ones(n2,1)/(2*alpha);
      ee = ones(n2,1)/alpha;
      count = 0; 
      for j = 1:n
         for i = 1:j
	    row = rr + (i-1)*n; col = cc + (j-1)*n;
	    if (i==j) 
               Fcell2{count+1} = spconvert([row,col,ee;n2,n2,0]); 
	       count = count+1; 
	    elseif (i < j)
               Fcell2{count+1} = spconvert([row,col,vv;n2,n2,0]); 
   	       count = count+1; 
            end
         end
      end
      idx = [1:n*(n+1)/2-1];
      Ftmp = svec(blk,Fcell(idx),1); Ftmp2 = svec(blk,Fcell2(idx),1); 
      Ft{1} = [Ft{1}, Ftmp{1}, Ftmp2{1}];
      bb = [bb; btmp(idx); btmp2(idx)];  
   end
%%******************************************************
