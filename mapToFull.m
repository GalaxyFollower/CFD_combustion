function [ Full ] = mapToFull( Dense,indexMap )

[y,x] = size(Dense);
l = max(indexMap);
if (y==x)
    Temp = zeros(l,x);
    Temp(indexMap(1:x),:)=Dense;
    Full = zeros(l,l);
    Full(:,indexMap(1:x))=Temp;
elseif y>x
        Full = zeros(l,x);
        Full(indexMap(1:y),:)=Dense;
elseif x>y
        Full = zeros(x,l);
        Full(:,indexMap(1:x))=Dense;
end

end

