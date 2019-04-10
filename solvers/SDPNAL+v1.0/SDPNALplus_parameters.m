%%**********************************************************************
%% SDPNALplus_parameters: parameters use in sdpnalplus
%%**********************************************************************
%% SDPNAL+: 
%% Copyright (c) 2014 by
%% Defeng Sun, Kim-Chuan Toh, Yancheng Yuan, Xinyuan Zhao
%% Corresponding author: Kim-Chuan Toh 
%% Last modified 2017/Oct/14
%%**********************************************************************

  function OPTIONS = SDPNALplus_parameters

  OPTIONS.tol        = 1e-6;  %accuracy tolerance for solving the SDP problem
  OPTIONS.maxiter    = 20000; %maximum number of outer iteration 
  OPTIONS.maxtime    = 10000; %maximum time allowed in seconds
  OPTIONS.tolADM     = 1e-4;  %%accuracry tolerance for using ADMM+ in Phase I of the algorithm
  OPTIONS.maxiterADM = 200;   %%maximum number of ADMM+ iterations to use in Phase I of the algorithm
  
  OPTIONS.stopoption = 1;     %(default) to terminate the solver when stagnation happens;
 %OPTIONS.stopoption = 0;     %otherwise.  
 
  OPTIONS.AATsolve.method = 'direct';    %(default) solve AAt*y=rhs by sparse Cholesky factor
 %OPTIONS.AATsolve.method = 'iterative'; %          solve AAt*y=rhs by diagonally preconditioned PSQMR
  OPTIONS.BBTsolve.method = 'iterative'; %(default) solve BBt*y2=rhs by diagonally preconditioned PSQMR
 %OPTIONS.BBTsolve.method = 'direct';    %          solve BBt*y2=rhs by sparse Cholesky factor
 
  OPTIONS.printlevel = 1;     %(default) level of details to print run history. 
                              % Value = 0 (less) , 1, 2 (more).
%%**********************************************************************

