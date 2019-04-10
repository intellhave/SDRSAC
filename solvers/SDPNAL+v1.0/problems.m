
   function [fname,fd] = problem
   
%%
%% theta: random
%%
   fname{11,1} = 'theta4';      
   fname{12,1} = 'theta42';    
   fname{13,1} = 'theta6';     
   fname{14,1} = 'theta62';    
   fname{15,1} = 'theta8';     
   fname{16,1} = 'theta82';    
   fname{17,1} = 'theta83';    

   fd(11:20) = 1*ones(1,10);  
%%
%% theta: Dimacs 
%%
   fname{21,1} = 'MANN-a27';      
   fname{22,1} = 'san200-0.7-1';  
   fname{23,1} = 'sanr200-0.7'; 
   fname{24,1} = 'c-fat200-1';    
   fname{25,1} = 'brock200-1';  
   fname{26,1} = 'brock200-4';  
   fname{27,1} = 'brock400-1';  
   fname{28,1} = 'keller4';     
   fname{29,1} = 'p-hat300-1';  

   fd(21:30) = 2*ones(1,10); 
%%
%% theta: 
%%
   fname{31,1} = '1dc.128';       
   fname{32,1} = '1et.128'; 
   fname{33,1} = '1tc.128'; 
   fname{34,1} = '1zc.128'; 
   fname{35,1} = '1dc.512';       
   fname{36,1} = '1et.512'; 
   fname{37,1} = '1tc.512'; 
   fname{38,1} = '1zc.512'; 
   
   fd(31:40) = 3*ones(1,10); 
%%
%% FAP: 
%%
   fname{41,1} = 'fap01'; 
   fname{42,1} = 'fap02'; 
   fname{43,1} = 'fap03'; 
   fname{44,1} = 'fap04'; 
   fname{45,1} = 'fap05'; 
   fname{46,1} = 'fap06'; 
   fname{47,1} = 'fap07'; 
   fname{48,1} = 'fap08'; 
   fname{49,1} = 'fap09'; 
   fname{50,1} = 'fap10'; 

   fd(41:50) = 4*ones(1,10); 
%%
%% QAP 
%%
   fname{51,1} = 'chr12a'; 
   fname{52,1} = 'chr12b';
   fname{53,1} = 'chr12c';
   fname{54,1} = 'had12';
   fname{55,1} = 'nug12';
   fname{56,1} = 'nug15';
   fname{57,1} = 'tai15a';
   fname{58,1} = 'tai15b';
   fname{59,1} = 'chr20a';
   fname{60,1} = 'chr20b';
   
   fd(51:60) = 5*ones(1,10); 
%%
%% BIQ
%%   
   fname{61,1} = 'be100.1'; 
   fname{62,1} = 'be120.3.1';    
   fname{63,1} = 'be150.3.1';    
   fname{64,1} = 'be200.3.1'; 
   fname{65,1} = 'be250.1'; 
   fname{66,1} = 'bqp100-1';  
   fname{67,1} = 'bqp250-1'; 
   fname{68,1} = 'bqp500-1'; 
   
   fd(61:70) = 6*ones(1,10);    
%%
%% sparse random SDPs of Rendl
%%
   fname{71,1} = 'Rn3m20p3'; 
   fname{72,1} = 'Rn3m25p3'; 
   fname{73,1} = 'Rn3m10p4'; 
   fname{74,1} = 'Rn4m30p3'; 
   fname{75,1} = 'Rn4m40p3'; 
   fname{76,1} = 'Rn4m15p4'; 
   fname{77,1} = 'Rn5m30p3'; 
   fname{78,1} = 'Rn5m40p3'; 
   fname{79,1} = 'Rn5m50p3'; 

   fd(71:80) = 7*ones(1,10);  
%%
%% NCM
%%
   fname{81,1} = 'NCM1n200H1'; 
   fname{82,1} = 'NCM1n200H2'; 
   fname{83,1} = 'NCM1n400H1'; 
   fname{84,1} = 'NCM1n400H2';
   fname{85,1} = 'NCM1n800H1';
   fname{86,1} = 'NCM1n800H2';
  
   fd(81:90) = 8*ones(1,10);       
%%
%% BIQ with additional valid inequality constraints
%%   
   fname{91,1} = 'be100.1'; 
   fname{92,1} = 'be120.3.1';    
   fname{93,1} = 'be150.3.1';    
   fname{94,1} = 'be200.3.1'; 
   fname{95,1} = 'be250.1'; 
   fname{96,1} = 'bqp100-1';  
   fname{97,1} = 'bqp250-1'; 
   fname{98,1} = 'bqp500-1'; 
   
   fd(91:100) = 9*ones(1,10);       
%%**************************************************************


