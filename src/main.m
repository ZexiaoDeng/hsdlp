clc ;
clear ;
close all ;


A = [1, 1, -2, 1; -1, -1, 2, -1; 2, -1, 4, 0; -1, 2, -4, 0];
b = [10, -10, 8, 4]';
c = [1, -2, 1, 0]';
c = -c;
[x, y, phi, z, w, psi, fval, status] = hsd(A, b, c);
x
f = -fval