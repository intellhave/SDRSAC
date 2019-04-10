%%****************************************
%%
%%
%%****************************************

  function [hh,mm,ss] = mytime(t);
 
  t = round(t); 
  h = floor(t/3600);
  m = floor(rem(t,3600)/60);
  s = rem(rem(t,60),60);
  
  hh = num2str(h);   
  if (h > 0) & (m < 10) 
     mm = ['0',num2str(m)];
  else
     mm = num2str(m);    
  end
  if (s < 10) 
     ss = ['0',num2str(s)];
  else
     ss = num2str(s);    
  end
%%**************************************** 

