%%***********************************************************************
%% rcpread: generate SDP data for Relxation of Clustering Problem
%%      
%% min { <-W,X> | X*e = e, <I,X> = K, X psd, X >=0}
%%
%%***********************************************************************
%% SDPNAL+ 
%% Copyright (c) 2014 by
%% Liuqin Yang, Defeng Sun, and Kim-Chuan Toh
%%***********************************************************************

   function [blk,At,C,b] = rcpread(fname,K)

   if (nargin == 1); K = 2; end; 
   datacell = importdata(fname);
   if ~iscell(datacell)
       if isfield(datacell,'data');    
           W = datacell.data;
       else
           W = datacell;
       end       
   else
       W = [];
       for i = 1:length(datacell)
           len = length(datacell{i});
           numstart = len+1;
           for j = 1:len
               if ~isnan(str2double(datacell{i}(j)))
                   numstart =  j;
                   break;
               end
           end
           numend = len;
           for j = numstart+1:len
               if isletter(datacell{i}(j))
                   numend =  j-1;
                   break;
               end
           end
           tmp = str2num(datacell{i}(numstart:numend));
           W = [W;tmp];
       end
   end
   [n,d] = size(W);
   blk{1,1} = 's'; blk{1,2} = n;   
   C{1} = -W*W';  
%%
   At{1} = [];
   for k = 1:n
      S = sparse([1:n],k,0.5,n,n,n);
      A = S+S';
      Asvec = svec(blk,A);
      At{1} = [At{1} Asvec];
   end
   A = sparse([1:n],[1:n],1,n,n,n);
   Asvec = svec(blk,A);
   At{1} = [At{1} Asvec];
   b = [ones(n,1);K];
%%******************************************************
