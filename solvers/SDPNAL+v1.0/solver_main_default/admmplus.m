%%**********************************************************************
%% A Convergent symmetric Gauss-Seidel based Proximal ADMM 
%% for solving SDP problems of the form:
%% (P) min { <C,X> | A(X)=b, X \in K and X \in P, B(X)-s=0, l<= s <=u}
%%     where K = psd cone, 
%%           P = { X | L <= X <=U } is a polyhedral set
%% (D) max  <b,y> - \delta_P^*(-Z2) -\delta_Q^*(-v)
%%     s.t. At(y) + Bt(y2) + Z1 + Z2 = C, y2-v = 0,
%%          Z1 positive semidefinite
%%**********************************************************************
%% SDPNAL+: 
%% Copyright (c) 2014 by
%% Defeng Sun, Kim-Chuan Toh, Yancheng Yuan, Xinyuan Zhao
%% Corresponding author: Kim-Chuan Toh 
%% Last modified 2016/06/01
%%**********************************************************************
   function [obj,X,s,y,Z1,Z2,y2,v,info,runhist] = ...
             admmplus(blk,AA,C,b,L,U,BB,l,u,OPTIONS,X,s,y,Z1,Z2,y2,v) 
            
   matlabversion = sscanf(version,'%f');
   matlabversion = matlabversion(1);   
   tmp = version('-release');
   matlabrelease = str2double(tmp(1:4));
   if (matlabversion < 8.3) && (matlabrelease < 2014)
      fprintf('\n Warning: SDPNAL+ is designed to run on Matlab 7.4 and above'); 
      fprintf('\n it may not run on older versions');   
   end
   randnstate = rng; 
   rng('default');    
   if (nargin < 10); OPTIONS = []; end
   if (nargin < 7); BB=[]; l=[]; u=[]; end
   if (nargin < 5); L=[]; U=[]; end
   if (nargin < 4); error(' number of input arguments must be at least 4'); end
   if (nargin>=7 && nargin<9); error(' BB,l,u not specified correctly'); end
   if (nargin>=5 && nargin<6); error(' L,U not specified correctly'); end   
%%-----------------------------------------------------------
%% parameters
%%-----------------------------------------------------------   
   stoptol  = 1e-6;
   sigma    = 1;  
   maxiter  = 20000; 
   maxtime  = 10000; 
   printlevel = 1;  
   rescale    = 1;   
   stopoption = 1;
   scale_data = 1; 
   AATsolve.method = 'direct'; 
   BBTsolve.method = 'iterative';
   beta = 1.618;
   if isfield(OPTIONS,'tol');      stoptol = OPTIONS.tol; end
   if isfield(OPTIONS,'sigma');    sigma = OPTIONS.sigma; end
   if isfield(OPTIONS,'maxiter');  maxiter = OPTIONS.maxiter; end 
   if isfield(OPTIONS,'maxtime');  maxtime = OPTIONS.maxtime; end    
   if isfield(OPTIONS,'printlevel'); printlevel = OPTIONS.printlevel; end
   if isfield(OPTIONS,'rescale');    rescale = OPTIONS.rescale; end    
   if isfield(OPTIONS,'stopoption'); stopoption = OPTIONS.stopoption; end   
   if isfield(OPTIONS,'AATsolve');   AATsolve = OPTIONS.AATsolve; end
   if isfield(OPTIONS,'BBTsolver');  BBTsolve = OPTIONS.BBTsolve; end
   if isfield(OPTIONS,'beta');       beta = OPTIONS.beta; end   
%%-----------------------------------------------------------
%% validate At, Bt, C
%%-----------------------------------------------------------
   Cinput = C; 
   [At,C,b,L,U,Bt,l,u,par] = validate(blk,AA,C,b,L,U,BB,l,u);
   m = length(b); m2 = par.m2; 
   blknum = size(blk,1);   
   Zero = ops(C,'zeros');
   par.existZ2 = (par.existL | par.existU);
%%
%% check initial variables
%% put SDP blocks before linear block, convert ublk to lblk
%%
   if ~exist('X','var') || isempty(X); X = Zero; end
   if ~exist('y','var') || isempty(y); y = zeros(m,1); end
   if ~exist('Z1','var') || isempty(Z1); Z1 = Zero; end
   if ~exist('s','var')  || isempty(s);  s = zeros(m2,1); end 
   if ~exist('y2','var') || isempty(y2); y2 = zeros(m2,1); end 
   if ~exist('Z2','var') || isempty(Z2); Z2 = Zero; end   
   if ~exist('v','var')  || isempty(v);  v = zeros(m2,1); end    
   for p=1:size(blk,1)
      pblk = blk(p,:); 
      if strcmp(pblk{1},'s') %% store in vectorized format
         if (size(X{p},1)==size(X{p},2)); X{p} = svec(pblk,X{p}); end
         if (size(Z1{p},1)==size(Z1{p},2)); Z1{p} = svec(pblk,Z1{p}); end
         if (size(Z2{p},1)==size(Z2{p},2)); Z2{p} = svec(pblk,Z2{p}); end
      end
   end         
   [blk,At,C,L,U,Bt,X,Z1,Z2,Zero,dim,sdpblkidx,ublkinfo] = ...
        convertdata(blk,At,C,L,U,Bt,X,Z1,Z2,Zero); 
   if (printlevel) && (norm(ublkinfo))
      fprintf('\n ##### convert ublk to lblk'); 
   end
   tstart = clock;  
%%-----------------------------------------------------------
%% scale data: Amap, ATmap, Bmap, BTmap
%%-----------------------------------------------------------
   Atorg = At; borg = b; Corg = C; Lorg = L; Uorg = U;
   Btorg = Bt; lorg = l; uorg = u;
   normCorg = 1+Fnorm(Corg); 
   normborg = 1+norm(borg);     
   if (m>0); existA = 1; else existA = 0; end
   if (m2>0); existB = 1; else existB = 0; end
   if (existA)
      normA2 = zeros(m,1); 
      for p=1:blknum          
          normA2 = normA2 + (sum(At{p}.*At{p}))';
      end
      normA = max(1,sqrt(normA2)); DA = spdiags(1./normA,0,m,m);
      for p = 1:blknum; At{p} = At{p}*DA; end
      b = b./normA; 
      if isfield(AA,'Amap');
         Amap = @(X) AA.Amap(X)./normA;
      else
         Amap = @(X) AXfun(blk,At,X); 
      end
      if isfield(AA,'ATmap');
         ATmap = @(y) AA.ATmap(y./normA);
      else
         ATmap = @(y) Atyfun(blk,At,y);
      end   
      y = y.*normA;
   else
      Amap = @(X) zeros(0,1); ATmap = @(y) Zero;
      normA = zeros(0,1); 
   end      
   if (existB)       
      normB2 = zeros(m2,1); 
      for p=1:blknum          
          normB2 = normB2 + (sum(Bt{p}.*Bt{p}))';
      end
      if (true)
         normB = 1; %%1: much better for ext-BIQ
         Bt = scaling(1./normB,Bt);          
      else
         normB = max(1,sqrt(normB2)); 
         DB = spdiags(1./normB,0,m2,m2);
         for p = 1:blknum; Bt{p} = Bt{p}*DB; end
      end
      if ~isempty(l); l = l./normB; end
      if ~isempty(u); u = u./normB; end
      if isfield(BB,'Bmap');
         Bmap = @(X) BB.Bmap(X)./normB;
      else
         Bmap = @(X) AXfun(blk,Bt,X); 
      end
      if isfield(BB,'BTmap');
         BTmap = @(y) BB.BTmap(y./normB);
      else
         BTmap = @(y) Atyfun(blk,Bt,y);
      end    
      y2 = y2.*normB; 
      v  = v.*normB; 
   else
      Bmap = @(X) zeros(0,1); BTmap = @(y) Zero;     
      normB = zeros(0,1);       
   end
%%
%% AATsolve
%%
   diagAAt = zeros(m,1);  
   for p = 1:blknum; diagAAt = diagAAt + sum(At{p}.*At{p})'; end
   par.diagAAt = max(1e-4,diagAAt); %%for SSN   
   if (existA)
      if isfield(AATsolve,'map'); %%important to check first     
         AATsolve.map = @(rhs) normA.*AATsolve.map(rhs.*normA); %%checked: correct
         AATsolve.method = 'map';         
      elseif strcmp(AATsolve.method,'direct')
         [LAAt,AAt] = cholAAt(blk,At,m);
         LAAt.isidentity = false;      
         dd = diag(AAt); 
         if (Fnorm(AAt-speye(m)) < 1e-12) 
            %%For BIQ, FAP and Theta(+) problems
            if (printlevel); fprintf('\n AAt is identity '); end
            LAAt.isidentity = true;
         elseif (Fnorm(AAt-spdiags(dd,0,m,m)) < 1e-15*norm(dd)) && printlevel
            if (printlevel); fprintf('\n AAt is diagonal '); end
         end
         AATsolve.direct = @(rhs) linsysolve(LAAt,rhs);
      elseif strcmp(AATsolve.method,'iterative')       
         AATsolve.psqmrmaxiter = 30;
         AATsolve.psqmrtolconst = 1e-3;
         AATsolve.precond = @(y) y./diagAAt; 
         AATsolve.matvec = @(y) matvecAAt(blk,At,y);         
      else
         error('AATsolve.method not defined'); 
      end
   else
      AATsolve.method = [];
   end   
%%
%% BBTsolve
%% 
   diagBBt = zeros(m2,1);    
   for p = 1:blknum; diagBBt = diagBBt + sum(Bt{p}.*Bt{p})'; end   
   vconst = 1;    
   if (existB)
      if strcmp(BBTsolve.method,'direct') 
         [LBBt,BBt] = cholAAt(blk,Bt,m2,vconst);
         LBBt.isidentity = false;      
         dd = diag(BBt); 
         if (Fnorm(BBt-speye(m2)) < 1e-12) 
            if (printlevel); fprintf('\n BBt is identity '); end
            LBBt.isidentity = true;
         elseif (Fnorm(BBt-spdiags(dd,0,m2,m2)) < 1e-15*norm(dd)) && printlevel
            if (printlevel); fprintf('\n BBt is diagonal '); end
         end
         BBTsolve.direct = @(rhs) linsysolve(LBBt,rhs); 
      elseif strcmp(BBTsolve.method,'iterative')
         opteigs.issym = 1; opteigs.maxit = 1000; opteigs.tol = 1e-12; 
         opteigs.v0 = randn(m2,1);      
         eignum = min(m2,5); %% TKC
         if (m2 > 1)
            [Q,Lam] = eigs(@(y) Bmap(BTmap(y)),m2,eignum,'lm',opteigs);   
            Lam = full(diag(Lam)); 
         else
            Q = 1; Lam = 1; 
         end
         vconst = sqrt(max(Lam))/4; 
         mu = Lam+vconst; idx = find(mu > 0);         
         mu = mu(idx); Q = Q(:,idx);          
         r2 = min(find(mu-mu(end) < 1e-8));
         muend = mu(r2);
         if (printlevel > 1) 
            fprintf('\n solve BBt by iterative solver: find the largest %2.1d eigenvalues',eignum);            
            fprintf(', use the largest %2.1d eigenvalues, maxeig=%3.2e,vconst=%5.4e',r2,max(Lam),vconst); 
         end         
         if (r2 > 1)
            Q1 = bsxfun(@times,Q(:,1:r2-1),sqrt(1/muend-1./mu(1:r2-1)).');
            BBTprecond = @(y) (1/muend)*y-Q1*(Q1'*y);  
         else
            BBTprecond = @(y) (1/muend)*y;  
         end
         BBTsolve.precond = BBTprecond; %%better       
         BBTsolve.matvec = @(y) matvecBBt(blk,Bt,y,vconst);         
      else
         error('BBTsolve.method not defined')         
      end      
   else
      BBTsolve.method = [];
   end  
   par.diagBBtpI  = max(1e-4,diagBBt+vconst*ones(m2,1)); %%for SSN; important to put here    
%%   
   if (scale_data==1)
      bscale = max(1,norm(b)); 
      Cscale = max(1,Fnorm(C));         
   elseif (scale_data==0)
      bscale = 1; 
      Cscale = 1;          
   end   
   if isfield(OPTIONS,'bscale'); bscale = OPTIONS.bscale; end 
   if isfield(OPTIONS,'Cscale'); Cscale = OPTIONS.Cscale; end 
   [b,l,u,L,U] = scaling(1/bscale,b,l,u,L,U);
   [C]   = scaling(1/Cscale,C);
   [X,s] = scaling(1/bscale,X,s./normB);
   [y,y2,v,Z1,Z2] = scaling(1/Cscale,normA.*y,normB.*y2,normB.*v,Z1,Z2);   
%%-----------------------------------------------------------
%% print header
%%-----------------------------------------------------------
   if (printlevel)
      fprintf('\n------------------------------------------------------------------------------------');      
      fprintf('\n ADMM+: A Convergent ADMM for SDP problems');           
      fprintf('\n------------------------------------------------------------------------------------');       
      if (printlevel)
         fprintf('\n  max dim. of sdp  var = %2.0d,',max(dim.sblk)); 
         fprintf('    num. of sdp  blk = %2.0d ',blknum);          
         if (m>0); fprintf('\n  num. of equality constraints = %2.0d',m); end   
         if (m2>0); fprintf('\n  num. of inequality constraints = %2.0d ',m2); end           
         if strcmp(AATsolve.method,'direct')
            fprintf('\n  solve AAt by direct solver'); 
         elseif strcmp(AATsolve.method,'iterative')
            fprintf('\n  solve AAt by iterative solver'); 
         elseif strcmp(AATsolve.method,'map')
            fprintf('\n  solve AAt by externally provided map');   
         else
            if (existA); error('  solver for AAt is not specified properly'); end
         end         
         if strcmp(BBTsolve.method,'direct')
            fprintf('\n  solve BBt by direct solver'); 
         elseif strcmp(BBTsolve.method,'iterative')
            fprintf('\n  solve BBt by iterative solver'); 
         elseif strcmp(BBTsolve.method,'map')
            fprintf('\n  solve BBt by externally provided map');   
         else
            if (existB); error('  solver for BBt is not specified properly'); end
         end                  
         fprintf('\n  scale_data = %2.0f',scale_data);
         fprintf('\n  bscale = %3.2e, Cscale = %3.2e',bscale,Cscale); 
         fprintf('\n------------------------------------------------------------------------------------');
      end
      fprintf('\n  iter| pinforg  dinforg  relgaporg|    <C,X>       <b,y>   |');
      fprintf(' time |  sigma  |');
      if (printlevel > 1)
          if (par.existZ2)
             fprintf(' beta| rankZ1,Z2| etaC1   etaC2   etaK1   etaK2 |');
          else
             fprintf(' beta|rankZ1| etaC1  etaK1 |');              
          end
          if (existB); fprintf(' CG |'); end
      end
      fprintf('\n')
   end
%%-----------------------------------------------------------    
%% main solver
%%-----------------------------------------------------------   
   par.normborg = normborg; par.normCorg = normCorg; 
   par.bscale = bscale; par.Cscale = Cscale;
   par.normA  = normA;    par.normB = normB;   
   par.existA = existA;   par.existB = existB;
   par.rescale = rescale; par.printlevel = printlevel;   
   par.stoptol = stoptol; par.stopoption = stopoption;
   par.vconst  = vconst;  
   par.maxiter = maxiter;
   par.maxtime = maxtime; 
   par.sigma   = sigma;
   par.beta    = beta;
   par.dim     = dim; par.m = m; par.m2 = m2;
   par.tstart  = tstart;
   [X,Y,s,shat,y,Z1,Z2,y2,v,info,runhist] = ...
      admmplus_main_default(blk,Amap,C,b,L,U,Bmap,l,u,...
      ATmap,BTmap,AATsolve,BBTsolve,At,Bt,par,X,s,y,Z1,Z2,y2,v);
%%-----------------------------------------------------------
%% recover orginal variables
%%-----------------------------------------------------------
   [X,Y,s,shat]   = scaling(bscale,X,Y,s.*normB,shat.*normB);
   [y,y2,v,Z1,Z2] = scaling(Cscale,y./normA,y2./normB,v./normB,Z1,Z2);
   [normX,normZ1,normZ2,normy,norms,normv] = Fnorm(X,Z1,Z2,y,s,v);
   [eigX,etaK1,etaC1,etaK2,etaC2,etaI] = ...
         competaorg(blk,Lorg,Uorg,lorg,uorg,X,Y,Z1,Z2,s,shat,v,1,1,1,1); 
   primobjorg = blktrace(blk,Corg,X); 
   dualobjorg = blktrace(blk,Y,Z2);
   if (existA)      
      dualobjorg = dualobjorg + full(borg'*y);   
   end
   if (existB)      
      dualobjorg = dualobjorg + full(s'*v);   
   end
   obj = [primobjorg,dualobjorg];
   Aty = ops(Atyfun(blk,Atorg,y),'+',ops(Z1,'+',Z2));   
   if (existB)
      Bty2 = Atyfun(blk,Btorg,y2);       
      BX = AXfun(blk,Btorg,X);
      etaN = norm(s-BX)/(1+norms);   
   else
      Bty2 = Zero; etaN = 0;
   end
   primfeasorg = max([norm(borg-AXfun(blk,Atorg,X))/normborg,etaN,etaI]);   
   dualfeasorg = max(Fnorm(ops(ops(Corg,'-',Aty),'-',Bty2))/normCorg,norm(v-y2)/(1+normv));
   ttime = etime(clock,tstart);      
   runhist.primfeasorg(info.iter) = primfeasorg;    
   runhist.dualfeasorg(info.iter) = dualfeasorg;   
   runhist.relgaporg(info.iter) = -diff(obj)/(1+sum(abs(obj)));    
   info.dim  = dim; info.m = m; info.m2 = m2;
   info.admtime   = ttime;   
   info.totaltime = ttime;   
   info.iterADM = info.iter;   
   info.primobjorg = primobjorg;  info.dualobjorg = dualobjorg;
   info.obj = [primobjorg,dualobjorg];   
   info.relgap = runhist.relgaporg(end);
   info.etaRp = primfeasorg; info.etaRd = dualfeasorg;   
   info.etaK1 = etaK1; info.etaK2 = etaK2; 
   info.etaC1 = etaC1; info.etaC2 = etaC2;   
   info.eta = max([primfeasorg,dualfeasorg,etaK1,etaK2,etaC1,etaC2]);
   info.mineigX = full(ops(eigX,'min')); info.maxeigX = full(ops(eigX,'max'));   
   info.PLX = min(0,full(ops(ops(X,'-',Lorg),'min')));     
   info.PUX = max(0,full(ops(ops(X,'-',Uorg),'max')));  
   if (existB)
      info.PlBX = min(0,min(BX-lorg));     
      info.PuBX = max(0,max(BX-uorg));  
   end
   info.normX  = normX; info.norms  = norms;
   info.normy  = normy; info.normZ1 = normZ1; info.normZ2 = normZ2;      
   if (max(primfeasorg,dualfeasorg) < 1.02*stoptol)   
      info.termcode = 0; 
      msg2 = sprintf(' SDP is solved to the required accuracy');
   elseif (primfeasorg < sqrt(stoptol) && dualfeasorg < 1.02*stoptol); 
      info.termcode = -1;
      msg2 = sprintf(' SDP is partially solved successfully: ');
      msg2 = [msg2, sprintf('primfeasorg < %3.2e, dualfeasorg < %3.2e',sqrt(stoptol),stoptol)];
   elseif (primfeasorg < 1.02*stoptol && dualfeasorg < sqrt(stoptol)); 
      info.termcode = -2;
      msg2 = sprintf(' SDP is partially solved successfully: ');
      msg2 = [msg2, sprintf('primfeasorg < %3.2e, dualfeasorg < %3.2e',stoptol,sqrt(stoptol))];
   elseif (primfeasorg < 10*stoptol && dualfeasorg < 10*stoptol); 
      info.termcode = -3;
      msg2 = sprintf(' SDP is partially solved successfully: ');
      msg2 = [msg2, sprintf('primfeasorg < %3.2e, dualfeasorg < %3.2e',10*stoptol,10*stoptol)];
   else
      info.termcode = 1;
      msg2 = sprintf(' SDP is not solved successfully');  
   end
   if (printlevel)  
      fprintf('\n--------------------------------------------------------------');
      fprintf('------------------');
   end
   info.msg2 = msg2;  
   if (printlevel) 
      if ~isempty(info.msg); fprintf('\n %s',info.msg); end
      if ~isempty(info.msg2); fprintf('\n %s',info.msg2); end    
      printinfo(info);
   end
   [X,Z1,Z2] = restoreback(blk,X,Z1,Z2,Cinput,ublkinfo,sdpblkidx);   
   rng(randnstate); 
%%**********************************************************************
%%**********************************************************************
%% convert lblk back to ublk; restore to matrix format
%%**********************************************************************   
   
   function [X,Z1,Z2] = restoreback(blk,X,Z1,Z2,Cinput,ublkinfo,sdpblkidx)

   if (norm(ublkinfo))
      for p = 1:length(ublkinfo)
         if (ublkinfo(p))
            len = blk{p,2}/2; 
            X{p} = X{p}(1:len)-X{p}(len+[1:len]); 
            Z1{p} = sparse(len,1); Z2{p} = sparse(len,1); 
         end
      end
   end
   blk(sdpblkidx,:) = blk;
   X(sdpblkidx) = X; Z1(sdpblkidx) = Z1; Z2(sdpblkidx) = Z2;  
   for p=1:size(blk,1)
      pblk = blk(p,:); 
      if strcmp(pblk{1},'s') && (size(Cinput{p},1)==size(Cinput{p},2))
         %% restore to matrix format
         X{p} = smat(pblk,X{p});
         Z1{p} = smat(pblk,Z1{p});
         Z2{p} = smat(pblk,Z2{p});                 
      end
   end         
%%***********************************************************************
%%***********************************************************************
