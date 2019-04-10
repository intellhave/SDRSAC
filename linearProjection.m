% Use LP to solve the Linear assignment problem
function [X] = linearProjection(Xp)
    
    Xp(Xp<0) = 0;
    N = size(Xp,1);
    n = round(sqrt(N));   
    
    f = -Xp;
    % Constraint matrix    
    A = [];
    b = [];
    
%     for i=1:N
%         ar = zeros(1,2*N);
%         ar(i) = 1; ar(N+i) = -1;
%         A = [A; ar]; b=[b; Xp(i)];
%         
%         ar(i) = -1; ar(N+i) = -1;
%         A = [A; ar]; b=[b; -Xp(i)];        
%     end
    
    % Now, set sum of rows and column to 1
    AA = []; bb=[];
    for i=1:n 
        ta = zeros(n,n); ta(i,:) =1;
        tar = ta(:)';
        AA = [AA; tar];  bb = [bb;1];
        AA = [AA; -tar]; bb = [bb;-1];
        
        ta = zeros(n,n); ta(:,i) =1;
        tar = ta(:)';
        AA = [AA; tar]; bb = [bb;1];
        AA = [AA; -tar]; bb = [bb;-1];
        
    end
    
    A=[A; AA];
    b=[b;bb];
    % Now, solve lp
   
    %[Xlp] = gurobiLinProg(f, A,b);
    [Xlp] = linprog(f, A,b, [], [], zeros(N,1), ones(N,1));
    
    X = reshape(Xlp(1:N),n,n);
    %disp(X);
    
    




end