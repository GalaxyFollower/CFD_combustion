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
t_end = 1;
delta_t = 10^(-5);
 
steps = round((t_end-t_0)/delta_t);
timeInterval = t_0:delta_t:t_end-delta_t;

delta_t=(t_end-t_0)/(steps-1);%0.001 second=1millisecond

% Number of degrees of freedom (number of nodes per length)
dimX=100;
dimY=100;
ndof=dimX*dimY;

%Number of eigenvectors:
N = 20;

%Variables saved every iteration for post-processing
maxfreq = []; % Saves which frequency was the maximum frequency
maxfreqampl = []; %Saves the amplitude of the maximum frequency
maxampl = []; % Saves the overall maximum amplitude of the solution
fftsave = []; % Saves the whole fft (within a certain spectrum)
sourcesave = []; %Saves the whole source for every run
sourceproduct = []; %For showing the Rayleigh Criterion

%% set up mesh
length_range = [ones(1,20)*0.336,linspace(0.336,0.700,40),linspace(0.700,0.336,40)];
tau=100; %* nr of delta_t
eta_deta = zeros(2*N,tau+1);%time delay determines how much initial conditions we need
%%
for l2 = length_range
    [X,Y,delta_X,delta_Y,C,nodeInfo,boundOrientation] = SetUpMesh(dimX,dimY,h1,h2,l1,l2,c1,c2,f,d);
    indexMap = getIndexMap(nodeInfo);
    
    disp('mesh is set up');
    
    % figure(1);
    % spy(nodeInfo);
    
    % Fill matrix A and vector B. Solve the linear system
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
    
    pos_flame=find(nodeInfo==66);%1D-nummering van Y(boven-onder) over X(links-rechts)
    pos_ref_point=find(nodeInfo==67);
    pos_measure_point=find(nodeInfo==68);
    
    %Initial condition:
    %eta_deta(1:N,tau+1) = 1; %op t=tau+1 exitatie 1 mode
    
    %Feedback constants
    n=10; %linear%
    k=34; %nonlinear%
    
    disp('Solving for following constant');
    l2
    %eta_deta = zeros(2*N,tau+1); %start from still initial condition
    %eta_deta(1,tau+1)=1;
    [eta_deta, full_source] = solveSystem(D,V,Grad_V_full,N,delta_t,steps,dB,source,nodeInfo,sourceTemplate,eta_deta,tau,n,k);
    
    P_punt=V(find(indexMap==pos_measure_point),:)*eta_deta(1:N,:);
    %fftstart = round(steps*4/5);
    fftstart = 1;
    pfft = abs(fft(P_punt(fftstart:end)));
    maxfreqampl = [maxfreqampl, max(pfft(1:10000))];
    maxfreq = [maxfreq,find(pfft(1:10000) == maxfreqampl(end))];
    maxampl = [maxampl, max(P_punt)];
    fftsave = [fftsave, pfft(1:10000)'];
    sourcesave = [sourcesave,full_source'];
    sourceproduct =  [sourceproduct, P_punt(tau+2:end)*cumtrapz(full_source)]; %Integrate because feedback is the source integrated.
    disp('time iteration done... starting extracting pressure');
    
    eta_deta(:,1:2*tau)=eta_deta(:,end-2*tau+1:end); %For hysteresis: make sure to start new simulation at with the old values.
    
end
%P=V*eta_deta(1:N,:);%(ndof*N)(N*steps)

%% Create Gif with sparse Matrix
%CreateSparseGif(P,X,Y,100,dimX,dimY,indexMap,nodeInfo);

%% Create plots of peaks
Fs = 1/(delta_t*(length(pfft))); %Needed for translating the fourier frequencies into actual frequencies
%FS is the amount of times the measured period occures within one unit of time

figure;
semilogy(length_range,maxfreqampl,'*');
title('amplitude of the maximum frequency in the fft');
ylabel('amplitude');
xlabel('l2');

figure;
plot(length_range(21:60),maxfreq(21:60)*Fs,'*b');
hold on;
plot(length_range(61:100),maxfreq(61:100)*Fs,'.r');
title('the dominating frequency in the fft');
ylabel('frequency');
xlabel('l2');

figure;
semilogy(length_range(21:60),maxampl(21:60),'*b');
hold on;
semilogy(length_range(61:100),maxampl(61:100),'.r');
title('global maximum of the solution');
ylabel('amplitude');
xlabel('l2');

figure;
maxfrequency = 2500;
surf(length_range(21:60),(1:maxfrequency)*Fs,log(fftsave(1:maxfrequency,(21:60))),'linestyle','none');
title('Plot showing the fft result for different parameters, length increasing');
xlabel('l2');
ylabel('fft frequency');
zlabel('magnitude');
colorbar;

figure;
maxfrequency = 2500;
surf(length_range(61:100),(1:maxfrequency)*Fs,log(fftsave(1:maxfrequency,(61:100))),'linestyle','none');
title('Plot showing the fft result for different parameters, length decreasing');
xlabel('l2');
ylabel('fft frequency');
zlabel('magnitude');
colorbar;

figure;
semilogy(length_range,delta_t*sourceproduct,'*');
title('Rayleigh Criterion');
ylabel('Integral value');
xlabel('l2');
%% time evolution of pressure at a point
P_punt=V(find(indexMap==pos_measure_point),:)*eta_deta(1:N,:);
%[pfft,f] = pwelch(P_punt(tau+80000:end),500,300,500,1/delta_t);%every point at time step of delta_t=10^-6 seconde

figure;
plot(timeInterval,P_punt(tau+2:end),'.');
xlabel('time [s]');
ylabel('Pressure[kg/m s^2]=[Pa]');
%
% figure(5);
% interest_range = 100;
% plot(f(1:interest_range),10*log10(pfft(1:interest_range)));
% xlabel('Frequency (Hz)');
% ylabel('Magnitude (dB)');
%
figure;
plot(eta_deta(1,:),eta_deta(N+1,:));
title('phase space of eigenvector one');
xlabel('function');
ylabel('derivative');
%
% figure;
% pf = fft(P_punt(end-20000:end));
% plot(abs(pf(1:100)))
%
figure;
dP_punt = diff(P_punt);
plot(P_punt(end-20000:end),dP_punt(end-20000:end));
title('Phase space of whole solution');
xlabel('Pressure');
ylabel('Time derivative of pressure');
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
% figure(7);
% title('Speed of sound over the domain');
% C(find(nodeInfo<0))=NaN;
% surf(X,Y,C);
%%
% timelength = delta_t*(steps*1/5);
% timescale = linspace(0,timelength,length(pfft));
% timescale = timescale(1:length(fftsave(:,1)));
