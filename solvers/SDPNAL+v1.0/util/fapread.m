%%*******************************************************
%% fapproblem:
%%
%% (primal prob.) min  Tr C*X
%%                s.t. diag(X) = e 
%%                    -Eij * X = 2/(k-1)  if (i,j) in  U
%%                    -Eij * X <= 2/(k-1) if (i,j) in GU
%%  where Eij = ei*ej' + ej*ei'
%%  b = [e; 2/(k-1)*ones(mU+mGU,1)]
%%  C = (1/2)*Diag(We) - (k-1)/(2k) *L(G) 
%%-------------------------------------------------------
%% [blk,Avec,C,b] = fapproblme(fname);
%%
%% SDPT3: version 3.0 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 2 Feb 01
%%*******************************************************

 function [blk,Avec,C,b,kparm] = fapread(fname)

%%
%% read fap data
%%
   if exist(fname)
      fid = fopen(fname,'r');
   elseif exist([fname,'.dat']); 
      fid = fopen([fname,'.dat'],'r');
   else 
      fprintf('** Problem not found. \n'); 
      blk = []; Avec = []; C = []; b = []; kparm = [];
      return;
   end
   [tmpr,count] = fscanf(fid,'%c');
   datavec = sscanf(tmpr,'%f'); clear tmpr;
   n = datavec(1); 
   numedges = datavec(2);
   kparm = datavec(3);  
   datavec = datavec(4:length(datavec)); 
   len = length(datavec);
   if (len ~= 3*numedges)   
      error(' fapread: numedges and data do not match.');   
   end
   I = datavec(1:3:len); 
   J = datavec(2:3:len); 
   w = datavec(3:3:len);    

   idxU = find(w==1000); 
   IU = I(idxU); JU = J(idxU); wU = w(idxU);
   U = spconvert([IU JU wU; n n 0]);
   U = U + U'; 
 
   idx2 = find(w~=1000); 
   I2 = I(idx2); J2 = J(idx2); w2 = w(idx2);       
   GU = spconvert([I2 J2 w2; n n 0]); 
   GU = GU + GU'; 
    
   fclose(fid);
%%
%% blk, Avec, C, b
%%    
    n = length(U); 
    mU = nnz(triu(U,1)); 
    mGU = nnz(triu(GU,1)); 
    m = mGU + mU + n; 

    b = [ones(n,1); (2/(kparm-1))*ones(mU+mGU,1)]; 
    LG = diag(GU*ones(n,1))-GU;
    C{1,1} = 0.5*(diag(GU*ones(n,1))) - ((kparm-1)/(2*kparm))*LG; 
    C{2,1} = zeros(mGU,1); 

    blk{1,1} = 's';  blk{1,2} = n;   
    blk{2,1} = 'l';  blk{2,2} = mGU; 
%%
%%
%%
    r2 = sqrt(2); 
    I = zeros(m,1); J = zeros(m,1); w = zeros(m,1); 
    cnt = 0; 
    e = [1:n]'; 
    I(1:n) = e.*(e+1)/2; 
    J(1:n) = e; 
    w = ones(n,1);  
    cnt = cnt+n;
    for i = 1:n 
        idx = find(U(i,i+1:n)); 
        idx = idx+i;      %% adjust index.  
        len = length(idx); 
        I(cnt+[1:len]) = i + idx.*(idx-1)/2; 
        J(cnt+[1:len]) = cnt+[1:len]'; 
        w(cnt+[1:len]) = -r2*ones(len,1); 
        cnt = cnt + len; 
    end  
    for i = 1:n
        idx = find(GU(i,i+1:n)); 
        idx = idx+i;      %% adjust index.  
        len = length(idx); 
        I(cnt+[1:len]) = i + idx.*(idx-1)/2; 
        J(cnt+[1:len]) = cnt+[1:len]'; 
        w(cnt+[1:len]) = -r2*ones(len,1); 
        cnt = cnt + len; 
    end
    Avec{1,1} = spconvert([I,J,w; n*(n+1)/2, m, 0]); 
%%
    I = [1:mGU]'; J = [m-mGU+1:m]';  w = ones(mGU,1); 
    Avec{2,1} = spconvert([I, J, w; mGU, m, 0]); 
%%
%%*******************************************************

