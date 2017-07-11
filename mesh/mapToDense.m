function [ Dense ] = mapToDense( Full,indexMap )

[y,x] = size(Full);
l = length(indexMap);

if (y==x)
    Temp = zeros(l,x);
    Dense = zeros(l,l);
    for i=1:l
        Temp(i,:)=Full(indexMap(i),:);
    end
    for i = 1:l
        Dense(:,i)=Temp(:,indexMap(i));
    end
elseif (x>y)
    Dense = zeros(y,l);
    for i=1:l
        Dense(:,i)=Full(:,indexMap(i));
    end
elseif (y>x)
        Dense = zeros(l,x);
    for i=1:l
        Dense(i,:)=Full(indexMap(i),:);
    end
end


end

