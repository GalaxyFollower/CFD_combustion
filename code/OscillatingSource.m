function [ sourcefun ] = OscillatingSource( Vn,N,dimY,dimX,indexMap,nodeInfo,omega )
sourceMat = zeros(2*N,1);

b=zeros(dimY*dimX,1);
b(find(nodeInfo==66))=1;
b = mapToDense(b,indexMap);
sourceMat(N+1:end) = Vn'*b;

sourcefun = @(Vn,Dn,N,t,eta,dimY,dimX,indexMap,nodeInfo) sourceMat*cos(omega*t);
end

