function [ eta_deta,total_source ] = solveSystem(D,V,Grad_V_full,N,delta_t,steps,dB,source,nodeInfo,sourceTemplate,eta_deta,tau,n,k)

theta = 0.5;


sysMat = zeros(2*N);
sysMat(1:N,N+1:end)=eye(N);
sysMat(N+1:end,1:N)=D;

A_star=eye(2*N)-delta_t*theta*sysMat;
rho=0.5;%[kg/m^3]
Q_ref=2000; %[2000Watt=2000 kg*m^2/s^3]
u_ref=10; %[m/s]

pos_ref_point=find(nodeInfo==67);
Grad_V_ref_point=Grad_V_full(pos_ref_point,:);
u=zeros(1,tau+1);%[m/s]
total_source = zeros(size(source));

for i=tau+1:steps+tau
    dB_i = dB + sourceTemplate*source(i-tau);
    feedback=-Q_ref/u_ref*n/rho*Grad_V_ref_point*eta_deta(1:N,i-tau)*(1-k*u(1,i-tau)^2);%(1x12)(12x1)
    dB_i = dB_i + feedback*sourceTemplate;
    total_source(i-tau) = source(i-tau)+feedback; %Save full source for analysis.
    b=zeros(2*N,1);
    b(N+1:end,1)=V'*dB_i;%one minus at the feedback is enough
    b_star=delta_t*b+eta_deta(:,i)+((1-theta)*sysMat*delta_t)*eta_deta(:,i);
    eta_deta(:,i+1)=A_star\b_star;
    u(1,i+1)=u(1,i)-delta_t/rho*Grad_V_ref_point*eta_deta(1:N,i);
end

end

