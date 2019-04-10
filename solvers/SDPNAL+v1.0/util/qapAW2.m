%%***********************************************************************
%% qapAW: generate SDP data for the problem QAP-AW1 in 
%%        the paper by Povh and Rendl.
%%
%%      [Y_{11} ... Y_{1n} x1]
%%  X = [  :    ...   :    : ]
%%      [Y_{n1} ... Y_{nn} xn]
%%      [ x1'   ...  xn'    1]
%%
%%***********************************************************************
%% SDPNAL+ 
%% Copyright (c) 2014 by
%% Liuqin Yang, Defeng Sun, and Kim-Chuan Toh
%%***********************************************************************

   function [blk,Ft,CC,bb] = qapAW(A,B,options)
   
   n = size(A,1);
   if (size(B,1)~= n); 
      error('A and B must have the same size'); 
   end
   if (nargin == 2); options = 2; end
   
   n2 = n*n + 1;
   blk{1,1} = 's'; blk{1,2} = n2;
%%
%% 
%%
   CC{1} = kron(B,A); 
   if norm(A-A','fro') | norm(B-B','fro')
      CC{1} = 0.5*(CC{1} + CC{1}'); 
   end
   CC{1} = [CC{1}, zeros(n*n,1); zeros(1,n*n), 0]; 
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
%%
%% add <ee^T,Y^{ij}> = 1 for j=1:n, i=1:j. 
%% 
   Fcell2 = cell(1,n*(n+1)/2);
   alpha = sqrt(n);
   btmp2 = ones(n*(n+1)/2,1)/alpha;      
   rr = [1:n]'*ones(1,n); rr = rr(:);  
   cc = ones(n,1)*[1:n];  cc = cc(:);
   vv = ones(n*n,1)/(2*alpha);
   ee = ones(n*n,1)/alpha;
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
%%
%% add diag(Y_{ii}) = xi for i=1:n 
%% Y_{n2,n2} = 1; 
%% 
   Fcell = cell(1,n2); 
   for k = 1:n*n
      Fcell{k} = spconvert([k,k,1; k,n2,-0.5; n2,k,-0.5; n2,n2,0]);
   end
   Fcell{n2} = spconvert([n2,n2,1]);
   Ftmp = svec(blk,Fcell,1); 
   Ft{1} = [Ft{1}, Ftmp{1}];
   bb = [bb; zeros(n*n,1); 1];
%%
%% add e'*xi = 1 for i= 1:n
%%   
   pp = [1:n]'; vv = 0.5*ones(n,1); 
   Fcell = cell(1,n); 
   for k = 1:n
      col = n2*ones(n,1); row = (k-1)*n + pp; 
      Fcell{k} = spconvert([row,col,vv; col,row,vv; n2,n2,0]);
   end
   Ftmp = svec(blk,Fcell,1); 
   Ft{1} = [Ft{1}, Ftmp{1}];
   bb = [bb; ones(n,1)]; 
%%******************************************************
