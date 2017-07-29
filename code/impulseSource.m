function [ sourceTerm ] = impulseSource( Vn,Dn,N,t,eta,dimY,dimX,indexMap,nodeInfo)

sourceMat = zeros(2*N,1);

b=zeros(dimY*dimX,1);
b(find(nodeInfo==66))=10;
b = mapToDense(b,indexMap);
sourceMat(N+1:end) = Vn'*b;

sourceTerm = sourceMat*exp(-t*20);

end

