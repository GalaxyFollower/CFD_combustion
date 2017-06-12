function [ indexMap ] = getIndexMap(nodeInfo )

[y,x] = size(nodeInfo);

indexMap = [];

for i = 1:x
    for j = 1:y
        nodeIndex = index(j,i,x,y);
        if (nodeInfo(j,i)>-1)
            indexMap = [indexMap;nodeIndex];
        end
    end
end
end

