function T = solveFVM(T_u, X, Y, boundary, TD, alpha, Tinf,Scheme,delta_t,dimX,dimY,theta)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File solveFVM.m
%
% This routine set up the linear system and solve it
%
% input
% T         Spatial Matrix T
% X         Matrix x coordinates 
% Y         Matrix y coordinates
% boundary  String vector. Boundary types.
% TD        Temperature for each boundary (if Dirichlet)
% alpha     convective heat transfer coefficient
% Tinf      Temperature of the surrouding fluid 
%
% output
% T         Temperature field

% Index maps the node position to the correct linear equation

index = @(ii, jj) ii + (jj-1) * dimY;


% B is the right-hand side of the linear system
B = zeros(dimY,dimX);


% set up the system matrix A
A = zeros(dimX * dimY);
for i = 1:dimY
    for j = 1:dimX
        % Fill the system matrix and the right-hand side for node (i,j)
        [output1,output2]=stamp(i, j, X, Y,TD,alpha,Tinf,boundary);
        A(index(i,j), :) = output1;%row vector
        B(i,j)=output2;%B has the grid-dimension
    end
end

%Make matrices sparse:
A = sparse(A);
B = sparse(B);

%%%%%%%%%%%%%%%%%CASE:Steady or unsteady%%%%%%%%%%%%%%%%%%%
% solve the linear system
if strcmp(boundary.west, 'Dirichlet') %BC at westside
    T_u(:,1)=TD.west;
end

if strcmp(boundary.east, 'Dirichlet')
    T_u(:,dimX)=TD.east;
end

if strcmp(boundary.north, 'Dirichlet')
    T_u(1,:)=TD.north;
end

if strcmp(boundary.south, 'Dirichlet')
    T_u(dimY,:)=TD.south;
end

switch Scheme
        
    case 'Steady' %linear to solve: Jacobi,GaussSeidel, SOR !don't use the 'iter' files, those are for given amount of iterations, not with tolerance
       save('A_matrix','A', 'B')
        tolerance = 0.0001;
        relaxation= 1.64;%fastest for SOR
        iterations=1000;
        Tinit=zeros(size(B(:)));
        T = Lin_solver(A,B(:),Tinit,'SOR',tolerance,relaxation,iterations); %B is a matrix, B(:) produces the columns in a vector
        T = reshape(T,dimY,dimX);
        
    case 'Unsteady Forward Euler' %explicit method, no need linear systems to solve
        %start iteration
        i=1;
        T_it=T_u(:);
        for t=0:delta_t:1
            T_it(:,i+1)=T_it(:,i)+delta_t*(A*T_it(:,i)-B(:));
            i=i+1;
        end  
        T=T_it;
        
    case 'Unsteady Theta'   %linear to solve: Jacobi,GaussSeidel, SOR
        A_star=sparse(eye(dimX*dimY)-delta_t*theta.*A);
        i=1;
        T_it=T_u(:);
        tolerance = 0.0001;
        relaxation= 1.64;%fastest for SOR
        iterations=10000;
        Tinit=zeros(size(B(:)));
        for t=0:delta_t:1%bigger delta_t when theta>0.5
            B_star=-B(:).*delta_t+(delta_t-theta*delta_t).*A*T_it(:,i)+T_it(:,i);
            T_it(:,i+1)= Lin_solver(A_star,B_star,Tinit,'JAC',tolerance,relaxation,iterations);
            i=i+1
        end
        T=T_it;
     
    case 'Unsteady RK4' %explicit method, no need linear systems to solve
        i=1;
        T_it=T_u(:);
        for t=0:delta_t:1
            T_p=T_it(:,i)+1/2*delta_t.*A*T_it(:,i)-B(:).*delta_t/2;
            T_pp=T_it(:,i)+1/2*delta_t.*A*T_p-B(:).*delta_t/2;
            T_ppp=T_it(:,i)+delta_t.*A*T_pp-B(:).*delta_t/2;
            T_it(:,i+1)=T_it(:,i)+1/6*delta_t.*(A*T_it(:,i)+2.*A*T_p+2.*A*T_pp+A*T_ppp)-B(:)*delta_t;
            i=i+1;
        end
        T=T_it;
end
        

