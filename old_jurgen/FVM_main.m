clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program to solve the 2D steady heat equation in a non-Cartesian Grid by
% the Finite Volumes Method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Authors: Deniz, Ranyu, Jurgen, wanqi, robin, ranyu

% Initialize variables

InitFVM
fprintf('shape=%s\n',shape);
fprintf('boundary.west=%s\n',boundary.west);
fprintf('boundary.south=%s\n',boundary.south);
fprintf('boundary.east=%s\n',boundary.east);
fprintf('boundary.north=%s\n',boundary.north);
fprintf('Scheme=%s\n',Scheme);

% set up the mesh

[X, Y] = setUpMesh(l, formfunction,dimX,dimY);

% Fill matrix A and vector B. Solve the linear system
T = solveFVM(T_u, X, Y, boundary, TD, alpha, Tinf,Scheme,delta_t,dimX,dimY,theta);

% Make some plots
VisTemperature(formfunction, X, Y, T, l,Scheme,dimX,dimY)








