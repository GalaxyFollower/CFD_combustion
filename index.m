function index = index(i,j,dimY)
%Index gives the place in the column vector where the value of the variable
%of the point in the grid is saved

index = i + (j-1) * dimY;

% if i==0
%     index = -1; %'North boundary'
% elseif i>dimY
%     index = -2; %'South boundary'
% elseif j == 0
%     index = -3; %'West boundary'
% elseif j > dimX
%     index = -4; %'East boundary'
% end


end

