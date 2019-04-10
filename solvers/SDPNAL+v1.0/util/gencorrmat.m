%%***********************************************************************
%% gencorrmat2: generate correlation matrix
%%              
%%  B = gencorrmat(n,options);
%% 
%%  options = 'Randcorr'
%%          = 'Randcorr2'
%%          = 'Diag'
%%          = 'AR1'
%%          = 'CompSym'
%%          = 'Rand'
%%***********************************************************************
%% SDPNAL+ 
%% Copyright (c) 2014 by
%% Liuqin Yang, Defeng Sun, and Kim-Chuan Toh
%%***********************************************************************

  function B = gencorrmat(n,options,alpha)

   rng('default')
   if nargin == 2; alpha = 0; end
  
   if strcmp(options,'Randcorr') || (options == 1)
      beta = 10^(-4/(n-1));
      xx = beta.^[0:n-1]; 
      xx = n*xx/sum(xx); 
      B = gallery('randcorr',xx);
   elseif strcmp(options,'Randcorr2') || (options == 2)
      n2 = n/2; 
      beta = 10^(-4/(n-1));
      xx = [beta.^[0:n2-1], zeros(1,n2)]; 
      xx = n*xx/sum(xx); 
      B = gallery('randcorr',xx);
   elseif strcmp(options,'Diag') || (options == 3)
      n2 = n/2; 
      tmp = [ones(n2), zeros(n2); zeros(n2), eye(n2)]; 
      B = diag(1e4*(2*rand(n,1)-1));
   elseif strcmp(options,'AR1') || (options == 4)
      tmp = -diag(ones(n,1)) - 0.8*diag(ones(n-1,1),-1); 
      Dinv = (1/0.01)*speye(n,n); 
      Binv = tmp'*Dinv*tmp; 
      B = inv(Binv); B = 0.5*(B+B'); 
      diagB = diag(B);
      B = diag(1./sqrt(diagB))*B*diag(1./sqrt(diagB)); 
   elseif strcmp(options,'CompSym') || (options == 5)
      sigma = 1; rho = 0.5; 
      dd  = sigma^2*(1-[0:n-1]*(rho^2)./(1+[0:n-1]*rho)); 
      tmp = zeros(n,n); 
      for j=2:n
         tmp(1:j-1,j) = rho./(1+[1:j-1]'*rho);
      end
      tmp = -eye(n) -tmp; 
      Dinv = diag(1./dd); 
      Binv = tmp'*Dinv*tmp;
      B = inv(Binv); B = 0.5*(B+B');    
      diagB = diag(B);
      B = diag(1./sqrt(diagB))*B*diag(1./sqrt(diagB)); 
   elseif strcmp(options,'Rand') || (options == 6)
      B = 2*rand(n)-1; 
      B = 0.5*(B+B');
   end
%% 
%% add noise
%%
    E = 2*(rand(n)-0.5); E = triu(E,1)+ triu(E,1)'; 
    B = (1-alpha)*B + alpha*E;
    B = min(B,1);
    B = max(-1,B);
    for i=1:n; B(i,i) = 1; end
%%****************************************************
