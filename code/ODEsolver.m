function [ t,eta ] = ODEsolver( Vn,Dn,eta0,N,tInterval,solverType,sourceFunction,dimX,dimY,indexMap,nodeInfo)

eta = zeros(2*N,1);

switch solverType
    case 'ode45' %Matlabs internal ode45 explicit solver
        sysMat = zeros(length(eta),length(eta));
        sysMat(1:N,N+1:end)=eye(N);
        sysMat(N+1:end,1:N)=Dn;
        
        fun = @(t,y) sysMat*y+sourceFunction(Vn,Dn,N,t,eta,dimY,dimX,indexMap,nodeInfo);
        
        [t,eta]  = ode45(fun,tInterval,eta0);
        
end


end

