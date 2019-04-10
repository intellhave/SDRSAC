function [A, b, C, X, y, Z] = rand_sdps(n, m, p, rand_seed)
% construct a primal dual pair of interior solutions to
%   P:  min <C,X>              D:  max b'y
%       st  A vec X = b              st  vec Z = vec C - A'y
%           X >= 0                       Z >= 0
% n ... primal dim
% m ... number of primal constraints
% A m by n*n; each row is vec( sym matrix)
% p ... size of nonzero submatrix forming A_i
% rand_seed ... random number initialisation
% call: [A, b, C, X, y, Z] = rand_sdps(n, m, p, rand_seed);

% initialize random number generator
rng(rand_seed);
 
% generate random Z and X, satisfying Z, X > pos. definite 
% Z
d = ones(n,1) + rand(n,1)*10;     
e = randn(n,ceil(n/4));
Z = diag( d) + e*e'/n;      
Z = (Z+Z')/2;

% X
d = ones(n,1) + rand(n,1)*10;
e = randn(n,ceil(n/4));
X = diag( d) + e*e'/n;      
X = (X+X')/2;

c_max = 100; % (the elements of Small_mat are between -5 and 5)
% generate random integer matrix A with sparse rows
A = sparse(n*n,m);
for j = 1:m
  if mod(j,10000)==0; fprintf(' j= %6.0d \n', j); end
% generate random selection of p rows/columns in perm
   perm = randperm(n); perm=perm(1:p); perm=sort(perm);
% generate a random matrix Small_mat of size p
   M = rand( p); M1 = triu( M,1);
%   Small_mat = M1+M1';
   Small_mat = diag(diag(M)) + M1 + M1';
% round so that numbers of Small_mat are between in (-c_max,....,cmax)
   Small_mat = round(-2*c_max * Small_mat + c_max);  
% remove main diag from Small_mat to reduce sparsity of A
   Small_mat = Small_mat - diag(diag( Small_mat));
%   Small_mat(1,p)=0; Small_mat(p,1)=0;  % more sparsity
% generate a matrix R with nonzero part in perm, given by Small_mat
   R = sparse(n,n); R(perm, perm) = Small_mat;
   R = sparse(R);
   A(:,j) = reshape(R, n*n,1);   % jth row of A
end
A=sparse(A);
A=A';

AAt = A*A';
perm = symamd( AAt);  
c = symbfact(AAt(perm,perm));
A=A(perm,:);

% select dual variables y normally distributed
% then b = A(X), and C = A'(y) + Z

y = randn(m,1);
b = A*reshape(X, n*n, 1);
C = reshape(A'*y, n, n) + Z;      
