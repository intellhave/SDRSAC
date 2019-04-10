%%**********************************************************************
%% SDPNAL+: A Majorized Semismooth Newton-CG
%%          Augmented Lagrangian Method to solve 
%%          primal and dual Semidefinite Programming
%%          of the form:
%% (P) min <C,X> 
%%     s.t. A(X) = b, B(X) = s, 
%%          X positive semidefinite
%%          X \in P:={X:L <= X <= U}, s \in Q:={s:l<= s<= u}
%% (D) max  <b,y> - \delta_P^*(-Z2) -\delta_Q^*(-v)
%%     s.t. At(y) + Bt(y2) + Z1 + Z2 = C, y2-v = 0,
%%          Z1 positive semidefinite
%%          
%% [obj,X,s,y,Z1,Z2,y2,v,info,runhist] = ...
%%      sdpnalplus(blk,AA,Cinput,b,L,U,BB,l,u,OPTIONS,X,s,y,Z1,Z2,y2,v)
%%
%% Input: blk,At,C,b (problem data whose input format is the same as in SDPT3) 
%%        L,U: follows the format of C. Set L or U = [] if no bound is imposed.
%%        Bt:  follows the format of At. Set Bt = [] if no inequality constraints exist.
%%        l,u: are vectors. Set l or u = [] if no bound is imposed. 
%%        OPTIONS: a structure array containing various parameters 
%%        Set OPTIONS = [] to use the default. For details, see SDPNALplus_parameters.m
%%
%% Output: obj  = [primal objective value, dual objective value]
%%         X    = a cell array containing the primal variable.
%%         y    = a vector of Lagrange multipliers to "A(X)=b"
%%         y2   = a vector of Lagrange multipliers corresponding to "B(X)=s"
%%         Z1,Z2= a cell array containing the dual slack variable.
%%         info = a structure array containing various information such as,
%%                iter    : number of outer iteration
%%                itersub : total number of sub-problems solved 
%%                obj     : objective value 
%%                relgap  : relative duality gap
%%                etaRp   : relative primal infeasibility 
%%                etaRd   : relative dual infeasibility 
%%                etaK1   : relative infeasibility for PSD-cone constraint
%%                etaK2   : relative infeasibility for bound constrints {X| L<= X<= U}
%%                etaC1   : relative complementarity for (X,Z1)
%%                etaC2   : relative complementarity for (X,Z2)
%%                eta     : max(etaRp,etaRd,etaK1,etaK2,etaC1,etaC2)
%%                PLX     : min(0,min(X-L))
%%                PUX     : max(0,max(X-U))
%%                termcode: termination code 
%%                          (0 : problem is solved to required tolerance, tol)
%%                          (-1: problem is partially solved to the required tolerance
%%                               with pinfeas < sqrt(tol) and dinfeas < tol)
%%                          ( 1: problem is not solved successful due to stagnation) 
%%                          ( 2: maximum number of iterations reached)
%%        runhist = a structure array containing the history of the run.
%%**********************************************************************
%% SDPNAL+: 
%% Copyright (c) 2014 by
%% Defeng Sun, Kim-Chuan Toh, Yancheng Yuan, Xinyuan Zhao
%% Corresponding author: Kim-Chuan Toh 
%% Last modified 2016/06/01
%%**********************************************************************
   warning off;
   clear all   
   restoredefaultpath;
   addpath(genpath(pwd)); 
   rng('default'); 
%%
   [fname,fd] = problems;     
   idxTheta1 = [11:17];       
   idxTheta2 = [21:29];
   idxTheta3 = [31:38];
   idxFAP    = [41:50];   
   idxQAP    = [51:60];
   idxBIQ    = [61:68];
   idxRand   = [71:79]; 
   idxNCM    = [81:86]; 
   idxextBIQ = [91:98];    
%%
   for k = [11]; [11,21,31,41,51,61,71,81,91];    
      probname = fname{k};   
      clear blk At C b L U Bt l u
      Bt = []; l = []; u = [];
      if (~isempty(probname))
         tstart=clock; 
         if (fd(k)==1)
            %% Theta+
            [blk,At,C,b] = read_sdpa([probname,'.dat-s']); 
            L = 0; U = [];
         elseif (fd(k)==2)
            %% Theta+
            [blk,At,C,b] = read_sedumi([probname,'.mat']); 
            L = 0; U = [];
         elseif (fd(k)==3)                
            %% Theta+
            [blk,At,C,b] = read_sdpa([probname,'.dat-s']); 
            L = 0; U = [];
         elseif (fd(k)==4)
            %% FAP
            [blk,At,C,b,L,U] = fapread_lu([probname]);                 
         elseif (fd(k)==5)
            %% QAP
            [AA,BB] = qapread([probname,'.dat']); 
            [blk,At,C,b,AAscale,BBscale] = qapAW(AA,BB,2);                
            L = 0; U = []; %% X >= 0
         elseif (fd(k)==6)
            %% BIQ
            [Q] = biqread([probname,'.sparse']); 
            [blk,At,C,b] = biq(Q); 
            L = 0; U = []; %% X >= 0               
         elseif (fd(k)==7)
            %% random sparse SDP
            n = str2num(probname(3))*100;
            idx1 = strfind(probname,'n');
            idx2 = strfind(probname,'m');
            idx3 = strfind(probname,'p');
            n = str2num(probname(idx1+1:idx2-1))*100;
            m = str2num(probname(idx2+1:idx3-1))*1000; 
            p = str2num(probname(idx3+1:end));             
            rand_seed = str2num([num2str(n),num2str(m/1000),num2str(p),'0']);
            [A,b,C] = rand_sdps(n,m,p,rand_seed); 
            K.s = n; [blk,At,C,b] = read_sedumi(A,b,C(:),K);
            L = []; U = []; %% no bound constraints on X                
         elseif (fd(k)==8)
            %% nearest (L1-error) correlation matrix problem
            idxH = find(probname=='H'); 
            idxN = find(probname=='n'); 
            n = sscanf(probname(idxN+1:idxH-1),'%d',idxH-idxN-1); 
            Hopt = sscanf(probname(end),'%d',1);
            matopt = sscanf(probname(4),'%d',1); 
            alpha = 0.05;
            B = gencorrmat(n,matopt,alpha); 
            if (Hopt == 1) 
               p = 0.01; 
            else
               p = 0.20; 
            end
            H0 = sprand(n,n,p); H0 = triu(H0) + triu(H0,1)';
            H1 = sprand(n,n,p); H1 = triu(H1) + triu(H1,1)';
            H2 = 100*H0 + 0.01*H1;
            H3 = rand(n,n); H3 = triu(H3) + triu(H3,1)';
            H3 = 0.1 + 10*H3;
            H = H2 + (1-spones(H2)).*H3;
            [blk,At,C,b] = NCM_Hnorm(B,H);
            L = []; U = []; %% no bound constraints on X      
         elseif (fd(k)==9)         
            [Q] = biqread([probname,'.sparse']); 
            [blk,At,A,C,b,Bt,Bstruct,l,u] = biq_addineq(Q,3);    
            probname = ['extended-',probname];  
            L = 0; U = [];
         end
         m = length(b); 
         ttime = etime(clock,tstart); 
         fprintf('\n\n %8s: ',probname); 
         fprintf(' read data time = %5.1f\n',ttime);
         OPTIONS.tol = 1e-6;     
         OPTIONS.stopoption = 0; 
         [obj,X,s,y,Z1,Z2,y2,v,info,runhist] = ...
             sdpnalplus(blk,At,C,b,L,U,Bt,l,u,OPTIONS);
         fprintf('\n press enter to continue\n');
         pause;             
      end 
   end 
%%******************************************************************
