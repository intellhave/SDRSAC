
[Introduction]

SDPNAL+ version 1.0: 
A MATLAB software for semidefinite programming with bound constraints
Copyright (c) 2014 by
Authors: Defeng Sun, Kim-Chuan Toh, Yancheng Yuan, Xinyuan Zhao
Corresponding author: Kim-Chuan Toh 
Past contributor: Liuqin Yang

This is software package for solving primal 
and dual SDP of the form: 
%%
%% (P) min { <C,X> | A(X)=b, X \in K and X \in P, B(X)-s=0, s\in Q}
%%     where K = psd cone, 
%%           P = { X | L <= X <= U} 
%%           Q = { s | l <= s <= u} 
%% (D) max  <b,y> - \delta_P^*(-Z2) -\delta_Q^*(-v)
%%     s.t. At(y) + Bt(y2) + Z1 + Z2 = C, y2-v = 0,
%%          Z1 positive semidefinite
%%
The main algorithm is a semi-smooth Newton CG augmented Lagrangian
method applied to (D); details can be found in the following 
reference: 

[Reference]
[1] L.Q. Yang, D.F. Sun, and K.C. Toh, 
SDPNAL+: a majorized semismooth Newton-CG augmented Lagrangian method 
for semidefinite programming with nonnegative constraints, 
Mathemtical Programming Computation, 7 (2015), pp. 331-366.
[2] X.Y. Zhao, D.F. Sun, and K.C. Toh, A Newton-CG augmented Lagrangian method 
for semidefinite programming, SIAM J. Optimization, 20 (2010), pp. 1737--1765. 

[Copyright] 
See Copyright.txt

--------------------------------------------------------------
[Installation] 
Run MATLAB in the directory SDPANL+v1.0.

After that, to see whether you have installed the software
correctly, type the following in MATLAB: 

>> SDPNALplus_Demo

or 
>> Example_NCM

or 
>> [blk,At,C,b] = read_sdpa('theta4.dat-s');
>> [obj,X,y,Z]  = sdpnalp(blk,At,C,b);
--------------------------------------------------------------