clear all
close all
clc
%% INIT
% Define dimension of the shape [m]
h1 = 0.200;
h2 =0.090;
l1=0.170;
l2=0.736;%max0.736,min 0.336
c1=0.080;
c2=0.180;
f=0.060; %flame length
d=0.080; %length between flame end and detection point

%time [s]
t_0 = 0;
t_end = 2;
delta_t = 10^(-5);

steps = round((t_end-t_0)/delta_t);
timeInterval = t_0:delta_t:t_end-delta_t;

tau=300; %* nr of delta_t
delta_t=(t_end-t_0)/(steps-1);%0.001 second=1millisecond

% Number of degrees of freedom (number of nodes per length)
dimX=100;
dimY=100;
ndof=dimX*dimY;

%Number of eigenvectors:
N = 2;

%Initial condition:
eta_deta = zeros(2*N,tau+1);%time delay determines how much initial conditions we need
%eta_deta(1:N,tau+1) = 1; %op t=tau+1 exitatie 1 mode

%Feedback constants
n=1000; %linear%
k=500; %nonlinear%

%% set up mesh
[X,Y,delta_X,delta_Y,C,nodeInfo,boundOrientation] = SetUpMesh(dimX,dimY,h1,h2,l1,l2,c1,c2,f,d);
indexMap = getIndexMap(nodeInfo);

disp('mesh is set up');

figure(1);
spy(nodeInfo);

%% Fill matrix A and vector B. Solve the linear system
% Solution of ddP-nabla.(c^2 nabla P)=q ==> ddP-AP=-B
%always at least timesteps 1/(2*f)=1/(2*omega/(2pi))
source_typ='WhiteNoise';

switch source_typ
    case 'cosinus'
        K=1;
        omega= sqrt(40);%omega=sqrt(-lambda) ; f=omega/(2pi)
        source=K*cos(omega*timeInterval);
    case 'WhiteNoise'
        fraction_nyq=1000/(steps/2);%we willen cut-off op 1000Hz
        source=WhiteNoise(steps,fraction_nyq,0,1);
    case 'nothing'
        source=zeros(1,steps);
end
sourceTemplate = zeros(ndof,1);
sourceTemplate(find(nodeInfo == 66))=1;
sourceTemplate = mapToDense(sourceTemplate,indexMap);

[A,B] = createSystemMatrices(dimX,dimY,X,Y,nodeInfo,boundOrientation,C);

disp('system matrices are assembled');

% Create dense matrices including only the nodes on the domain.
dA = sparse(mapToDense(A,indexMap));
dB = sparse(mapToDense(B(:),indexMap));

% Clear unneeded variables (just for memory)
clear A;

%Eigenvektors
[V,D] = eigs(dA,N,-0.01);

Grad_V_full=[];
for i=1:N
    V_i=V(:,i);
    V_full_i=mapToFull(V_i,indexMap,ndof);
    V_full_i=reshape(V_full_i,[dimY,dimX]);
    [dV_dx,dV_dy]=gradient(V_full_i,delta_X,delta_Y);
    dV_dx(find(nodeInfo==1))=0;
    dV_dx(find(nodeInfo==5))=0;
    dV_dx(find(nodeInfo==7))=0;
    dV_dx(find(nodeInfo==9))=0;
    dV_dx(find(nodeInfo==10))=0;
    dV_dx(find(nodeInfo==11))=0;
    dV_dx(find(nodeInfo==12))=0;
    dV_dx_d=mapToDense(dV_dx(:),indexMap);
    Grad_V_full=[Grad_V_full dV_dx(:)];
end

disp('Eigendecomposition of system is done...');

[eta_deta, full_source] = solveSystem(D,V,Grad_V_full,N,delta_t,steps,dB,source,nodeInfo,sourceTemplate,eta_deta,tau,n,k);

pos_flame=find(nodeInfo==66);%1D-nummering van Y(boven-onder) over X(links-rechts)
pos_ref_point=find(nodeInfo==67);
pos_measure_point=find(nodeInfo==68);

disp('time iteration done... starting extracting pressure');

%P=V*eta_deta(1:N,:);%(ndof*N)(N*steps)

%% Create Gif with sparse Matrix
%CreateSparseGif(P,X,Y,100,dimX,dimY,indexMap,nodeInfo);

%% time evolution of pressure at a point
P_punt=V(find(indexMap==pos_measure_point),:)*eta_deta(1:N,:);
%[pfft,f] = pwelch(P_punt(tau+80000:end),500,300,500,1/delta_t);%every point at time step of delta_t=10^-6 seconde

figure(4);
plot(timeInterval,P_punt(tau+2:end),'.');
xlabel('time [s]');
ylabel('Pressure[kg/m s^2]=[Pa]');

figure;
plot(eta_deta(1,:),eta_deta(N+1,:));
title('phase space of eigenvector one');
xlabel('function');
ylabel('derivative');
3

figure;
pf = fft(P_punt(end-2000:end));
plot(abs(pf(1:100)))
4

figure;
dP_punt = diff(P_punt);
plot(P_punt(end-2000:end),dP_punt(end-2000:end));
title('Phase space of whole solution');
%% make some plots of eigenvectors A
% [Vn,Dn] = eigs(dA,30,-0.01);
% figure(6);
% for i = 1:6
%     VFull = mapToFull(Vn(:,i),indexMap,ndof);
%     VFull = reshape(VFull,[dimY,dimX]);
%     VFull(find(nodeInfo<0))=NaN;
%     subplot(3,2,i);
%     surf(X,Y,VFull);
%     title(Dn(i,i));
% end
%     
figure(7);
title('Speed of sound over the domain');
C(find(nodeInfo<0))=NaN;
surf(X,Y,C);
xlabel('x');
ylabel('y');
zlabel('C');