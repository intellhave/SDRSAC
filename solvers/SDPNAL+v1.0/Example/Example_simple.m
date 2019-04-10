%%********************************************************************    
%% This is an example for the interface
%% This example is to solve the following SDP
%%
%% min trace(X1)+trace(X2)+sum(X3)
%% s.t. -X1(1,2)+2X2(3,3)+2X3(2)=4
%%       2X1(2,3)+X2(4,2)-X3(4) = 3
%%       2<= -X1(1,2)-2X2(3,3)+2X3(2)<=7
%%       X1 \in S^{6}_{+}, X2 \in R^{5*5}, X3 \in R^{7}_{+}
%%       0<=X1<=10, 0<=X2<=8
%% SDPNAL+: 
%% Copyright (c) 2017 by
%% Yancheng Yuan, Kim-Chuan Toh, Defeng Sun and Xinyuan Zhao
%%********************************************************************    
clear all;
n1 = 6; n2 = 5; n3 = 7;
%% Initialize model
model = ccp_model('Example_simple');
    X1 = var_sdp(n1,n1);
    X2 = var_nn(n2, n2);
    X3 = var_nn(n3);
    model.add_variable(X1, X2, X3);
    model.minimize(trace(X1)+trace(X2)+sum(X3));
    options=3;    
    if (options==1)
       model.add_affine_constraint(-X1(1,2)+2*X2(3,3)+2*X3(2)==4);
       model.add_affine_constraint(2*X1(2,3)+X2(4,2)-X3(4) == 3);
    else
       %%two other ways to add the above equality constraints
       A1 = {sparse(n1,n1), sparse(n1,n1)}; 
       A2 = {sparse(n2,n2), sparse(n2,n2)}; A3 = sparse(2,n3);
       A1{1}(1,2)= -1;  A1{2}(2,3)= 2; 
       A2{1}(3,3)= 2;   A2{2}(4,2)= 1; 
       A3(1,2) = 2;     A3(2,4) = -1;        
       b = [4;3]; 
       if (options==2)
          model.add_affine_constraint(A1*X1+A2*X2+A3*X3 == b);       
       elseif (options==3)
          blk{1,1} = 's'; blk{1,2} = n1;
          %A1vec = [svec(blk,A1{1}),svec(blk,A1{2})]';
          A1vec = [A1{1}(:), A1{2}(:)]'; 
          %model.add_affine_constraint(A1vec*map_svec(X1)+A2*X2+A3*X3 == b);  
          model.add_affine_constraint(A1vec*map_vec(X1)+A2*X2+A3*X3 == b);
       end
    end
    model.add_affine_constraint(2<=-X1(1,2)-2*X2(3,3)+2*X3(2)<=7);    
    model.add_affine_constraint(0<= X1 <=10);
    model.add_affine_constraint(X2<=8);
    model.add_affine_constraint(map_diag(X1)==ones(n1,1));
    model.setparameter('BBTsolve.method','direct');
model.solve;
%%*********************************************************************