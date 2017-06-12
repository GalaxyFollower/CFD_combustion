clear all;
close all

%% INIT
% Define dimension of the shape
h = 10;
l1=5;
l2=10;
c1=2;
c2=5;

% Number of degrees of freedom (number of nodes per length)
dimX=60;
dimY=60;

% % Parameter for Robin BC
% alpha = 2;
% Tinf = 90;

% Boundary conditions 
    %0 == dirichlet with in second column the pressure (in bar)
    %1 == Neumann with in second column the pressure flux (in bar)
    %every row= boundary type: row1=west boundary         row9=Northwestcorner
    %                          row2=east boundary         row10=Northeastcorner
    %                          row3=north boundary        row11=Southwestcorner
    %                          row4=south boundary        row12=Southeastcorner
    %                          row5=east_inner boundary   row13=INW
    %                          row6=north_inner boundary  row14=INE
    %                          row7=west_inner boundary   row15=ISW
    %                          row8=south_inner boundary  row16=ISE
    
boundParameters = [0 100;1 0;1 0;1 0;1 0;1 0;1 0;1 0;1 0;1 0;1 0;1 0;1 0;1 0;1 0;1 0];

%%Scheme...Steady for now

%% Main

%set up mesh
[X,Y,nodeInfo,boundOrientation] = SetUpMesh(dimX,dimY,h,l1,l2,c1,c2);
indexMap = getIndexMap(nodeInfo);

figure(1);
spy(nodeInfo);

% Fill matrix A and vector B. Solve the linear system
[A,B] = createSystemMatrices(dimX,dimY,X,Y,nodeInfo,boundOrientation,boundParameters);
figure(2);
spy(A);
figure(3);
spy(B);

% Create dense matrices including only the nodes on the domain.
dA = mapToDense(A,indexMap);
dB = mapToDense(B(:),indexMap);

% Create eigenvectors V and eigenvalues D
[V,D] = eigs(dA,10,'sm');

    
%make some plots
...

for i = 1:6
    VFull = mapToFull(V(:,i),indexMap);
    VFull = reshape(VFull,[dimY,dimX]);
    VFull(find(nodeInfo<0))=NaN;
    subplot(3,2,i);
    surf(X,Y,VFull);
end

P = V*inv(D)*V'*dB;

PFull = mapToFull(P,indexMap);
PFull = reshape(PFull,[dimY,dimX]);
PFull(find(nodeInfo<0))=NaN;
figure;
surf(X,Y,PFull);
figure;
contour(X,Y,PFull);
