%%*********************************************************************
%% Solve FAP in using interface(using loops)
%% to solve the following FAP problem
%% max <((k-1)/2k)*L(G,W)-1/2*Diag(We),X>
%% s.t. diag(X)==1 ; X positive semidefinite
%%      <-E^{i,j},X>==2/(k-1)  for (i,j) in U_set
%%      <-E^{i,j},X> <= 2/(k-1)  for (i,j) in Edge_set\U_set
%% SDPNAL+: 
%% Copyright (c) 2017 by
%% Yancheng Yuan , Kim-Chuan Toh, Defeng Sun and Xinyuan Zhao
%%*********************************************************************

%% read FAP data: U_set, Edge_set
 clear all;
 fname = 'fap08';
 if exist(fname)
    fid = fopen(fname,'r');
 elseif exist([fname,'.dat']); 
    fid = fopen([fname,'.dat'],'r');
 else 
    error('** Problem not found. \n'); 
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
 idxE = find(w~=1000); 
 IE = I(idxE); JE = J(idxE); wE = w(idxE);
 fclose(fid);
 
 %% model the problem using the interface
 GE = spconvert([IE JE wE; n n 0]); 
 GE = GE + GE';
 GU = spconvert([IU JU wU; n n 0]);
 GU = GU + GU';
 LG = diag(GE*ones(n,1))-GE;
 C =  -0.5*diag(GE*ones(n,1)) + ((kparm-1)/(2*kparm))*LG;  
 
 model = ccp_model('Example_FAP');
 X = var_sdp(n, n);
 model.add_variable(X);
 model.maximize(inprod(C,X));
 model.add_affine_constraint(map_diag(X) == ones(n,1));
 const = 2/(kparm-1);
 for k=1:length(IU)
     Eij = spconvert([IU(k),JU(k),1; JU(k),IU(k),1; n,n,0]);
     model.add_affine_constraint(inprod(-Eij, X) == const);
 end
 for k=1:length(IE)
     Eij = spconvert([IE(k),JE(k),1; JE(k),IE(k),1; n,n,0]);
     model.add_affine_constraint(inprod(-Eij, X) <= const);
 end
 model.solve;
 %%*********************************************************************