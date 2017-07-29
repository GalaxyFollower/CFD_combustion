function [ X,Y,delta_X,delta_Y,C,nodeInfo,boundOrientation ] = SetUpMesh(dimX,dimY,h1,h2,l1,l2,c1,c2,f,d)
%   Set up outer grid
%   dimX, dimY: The amount of grid points in the directions X and Y.
%   h,l1,l2,c1,c2: Geometric parameters.
%
%   X,Y fields needed for plotting
%   C:speed of sound grid
%   nodeInfo: Holds information about each node. positive integers refer to boundaries, 0 is an inner node, -1 is an outer node, boundary nodes:
%   see grid, 66=flame, 67=ref_point
%   boundOrientation: Holds information about the boundary's orientation
%                      (1:W,2:E,3:N,4:S,5:NW,6:NE,7:SW,8:SE) !!Innercorners now as sides, otherwise(9:INW,10:INE,11:ISW,12:ISE)
%%%%%%%%%%%%%%%%%%%%%%%%%%GEOMETRY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      -l1-----            ----l2-------
%                  -c2-
%  |  ( 9)333(10)      
%  |  1         5      ( 9)733333333(10)   |     |
%  |  1         5      7 (66)          2   c1    |
%  |  1        (5)6666(7)(66)          2   |     |
%  h1 1               (67)(66)(68)     2         h2
%  |  1        (5)8888(7)(66)          2         |  
%  |  1         5      7 (66)          2         |
%  |  1         5      (11)444444444(12)         |
%  |  (11)444(12)      


length_X = l1+l2+c2;
length_Y = h1;

delta_X = length_X/(dimX-1);
delta_Y = length_Y/(dimY-1);

x = linspace(0,length_X,dimX);
y = linspace(0,length_Y,dimY);

[X,Y] = meshgrid(x,y);

%% Build nodeInfo
nodeInfo = zeros(dimY,dimX);
indexLeft = ceil(l1/(delta_X));%net linkse boundary
indexRight = ceil((l1+c2)/delta_X);%net rechtse boundary
indexFlameEnd=ceil((l1+c2+f)/delta_X);
indexDetection=ceil((l1+c2+f+d)/delta_X);

indexBottomStart = ceil((h1-c1)/delta_Y);%net boven onderste boundary
indexTopEnd = ceil(c1/delta_Y);%net bovenste boundary

indexBottomStart_2=ceil((h1-(h1-h2)/2)/delta_Y);%net boven onderste boundary
indexTopEnd_2=ceil((h1-h2)/(2*delta_Y));%net bovenste boundary

%Flame+ref_point
nodeInfo(indexTopEnd_2+1:indexBottomStart_2,indexRight+1:indexFlameEnd)=66;
nodeInfo(ceil((indexTopEnd_2+indexBottomStart_2)/2),indexRight-5)=67;%feedback point, at start of combustion chamber
nodeInfo(ceil((indexTopEnd_2+indexBottomStart_2)/2),indexDetection)=68;

%Boundaries
nodeInfo(1,1:indexLeft)=3;
nodeInfo(indexTopEnd_2,indexRight:end)=3;

nodeInfo(end,1:indexLeft)=4;
nodeInfo(indexBottomStart_2+1,indexRight:end)=4;

nodeInfo(:,1)=1;
nodeInfo(indexTopEnd_2:indexBottomStart_2,end)=2;

% Make mesh indentations
tempblock1 = -1*ones(indexTopEnd,indexRight-indexLeft+1);
tempblock2= -1*ones(indexTopEnd_2-1,dimX-indexRight+1);

topblock = tempblock1;
topblock(end,:) = 6;
topblock(:,1)=5;%so standing sides after bottom
topblock(:,end)=7;

bottomblock = tempblock1;
bottomblock(1,:) = 8;%so standing sides after top
bottomblock(:,1)=5;
bottomblock(:,end)=7;

nodeInfo(1:indexTopEnd,indexLeft:indexRight)=topblock;
nodeInfo(indexBottomStart+1:end,indexLeft:indexRight)=bottomblock;
nodeInfo(1:indexTopEnd_2-1,indexRight:end)=tempblock2;
nodeInfo(indexBottomStart_2+2:end,indexRight:end)=tempblock2;

nodeInfo(1,1)=9;%NW corner
nodeInfo(indexTopEnd_2,indexRight)=9;%NW corner
%nodeInfo(indexTopEnd,indexRight)=13;%INW corner

nodeInfo(1,indexLeft)=10;%NE corner
nodeInfo(indexTopEnd_2,end)=10;%NE corner
%nodeInfo(indexTopEnd,indexLeft)=14;%INE corner

nodeInfo(end,1)=11;%SW corner
nodeInfo(indexBottomStart_2+1,indexRight)=11;%SW corner
%nodeInfo(indexBottomStart+1,indexRight)=15;%ISW corner

nodeInfo(end,indexLeft)=12;%SE corner
nodeInfo(indexBottomStart_2+1,end)=12;%SE corner
%nodeInfo(indexBottomStart+1,indexLeft)=16;%ISE corner

%% boundary orientation
boundOrientation=[1 1;2 2;3 3;4 4;5 2;6 3;7 1;8 4;9 5;10 6;11 7;12 8]; %[13 9;14 10;15 11;16 12]; Inner corners as sides

%% Speed of Sound
C=300*(2+tanh((X-(l1+c2))));%[mm/ms]
%C=300.*ones(size(X));
end