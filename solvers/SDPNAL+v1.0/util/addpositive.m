%%***********************************************************************
%% addpositive: add the constraints that X_ij is nonnegative. 
%%
%% [blk2,At2,C2,b2] = addpositive(blk,At,C,b);    
%%
%%***********************************************************************
%% SDPNAL+ 
%% Copyright (c) 2014 by
%% Liuqin Yang, Defeng Sun, and Kim-Chuan Toh
%%***********************************************************************

   function [blk2,At2,C2,b2] = addpositive(blk,At,C,b)

   m = length(b);
   alpha = 0.1;  
   
   dim = zeros(1,2); numblk = zeros(1,2); 
   sdpblkidx = [];    
   for p = 1:size(blk,1)
      pblk = blk(p,:); 
      if strcmp(pblk{1},'s')
         dim(1) = dim(1) + sum(pblk{2});
         numblk(1) = numblk(1) + length(pblk{2}); 
         sdpblkidx = [sdpblkidx, p]; 
      elseif strcmp(pblk{1},'l')
         dim(2) = dim(2) + sum(pblk{2});
         numblk(2) = numblk(2) + length(pblk{2}); 
      end
   end   
   sdpblkidx = [sdpblkidx, setdiff([1:size(blk,1)],sdpblkidx)];
   blk2 = blk(sdpblkidx,:); At2 = At(sdpblkidx); C2 = C(sdpblkidx); 
%%
   npos = 0; 
   numsblk = 0;   
   for p = 1:size(blk2,1)      
      pblk = blk2(p,:);
      if strcmp(pblk{1},'s')
         npos = npos + sum(pblk{2}.*(pblk{2}+1))/2;
      end
   end
%%
   len = 0; 
   for p = 1:size(blk2,1)      
      pblk = blk2(p,:);
      n = sum(pblk{2});   
      if strcmp(pblk{1},'s')
         n2 = sum(pblk{2}.*(pblk{2}+1))/2; 
         tmp = [sparse(n2,len),alpha*speye(n2),sparse(n2,npos-len-n2)]; 
         At2{p,1} = [At2{p,1}, tmp];
         len = len + n2; 
      else         
         At2{p,1} = [At2{p,1}, sparse(n,npos)];
      end
   end
   numblk = size(blk2,1); 
   blk2{numblk+1,1} = 'l'; blk2{numblk+1,2} = npos; 
   At2{numblk+1,1} = [sparse(npos,m), -alpha*speye(npos)];
   b2 = [b; zeros(npos,1)];
   C2{numblk+1,1} = zeros(npos,1);
%%***********************************************************
