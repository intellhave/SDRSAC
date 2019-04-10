function [R, t] = read4PCSResults(filename)

 M = dlmread(filename);
 R = M(1:3, 1:3);
 t = M(1:3, 4);
 



end