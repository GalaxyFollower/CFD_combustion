function [ X,Y,nodeInfo,boundOrientation ] = SetUpMesh(dimX,dimY,h,l1,l2,c1,c2)
%   Set up outer grid
%   dimX, dimY: The amount of grid points in the directions X and Y.
%   h,l1,l2,c1,c2: Geometric parameters.
%
%   X,Y fields needed for plotting
%   nodeInfo: Holds information about each node. positive integers refer to
%   boundaries, 0 is an inner node, -1 is an outer node.
%   boundOrientation: Holds information about the boundary's orientation
%   (1:W,2:E,3:N,4:S,5:NW,6:NE,7:SW,8:SE,9:INW,10:INE,11:ISW,12:ISE)
%%%%%%%%%%%%%%%%%%%%%%%%%%GEOMETRY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   -l1-------      ----l2----------
%              -c2-
%   ( 9)333(10)    ( 9)733333333(10)   |     |
%   1         5    7               2   c1    |
%   1      (14)6666(13)            2   |     |
%   1                              2         h
%   1      (16)8888(15)            2         |  
%   1         5    7               2         |
%   (11)444(12)    (11)444444444(12)         |


length_X = l1+l2+c2;
length_Y = h;

delta_X = length_X/(dimX-1);
delta_Y = length_Y/(dimY-1);

x = linspace(0,length_X,dimX);
y = linspace(0,length_Y,dimY);

[X,Y] = meshgrid(x,y);

%Build nodeInfo
nodeInfo = zeros(dimY,dimX);
nodeInfo(1,:)=3;
nodeInfo(end,:)=4;
nodeInfo(:,1)=1;
nodeInfo(:,end)=2;

%Make mesh indentations
indexLeft = ceil(l1/(delta_X));
indexRight = ceil((l1+c2)/delta_X);
indexBottomStart = ceil((h-c1)/delta_Y);
indexTopEnd = ceil(c1/delta_Y);

tempblock = -1*ones(indexTopEnd,indexRight-indexLeft+1);
tempblock(:,1)=5;
tempblock(:,end)=7;

topblock = tempblock;
topblock(end,:) = 6;

bottomblock = tempblock;
bottomblock(1,:) = 8;

nodeInfo(1:indexTopEnd,indexLeft:indexRight)=topblock;
nodeInfo(indexBottomStart+1:end,indexLeft:indexRight)=bottomblock;

nodeInfo(1,1)=9;%NW corner
nodeInfo(1,indexRight)=9;%NW corner
nodeInfo(indexTopEnd,indexRight)=13;%INW corner

nodeInfo(1,indexLeft)=10;%NE corner
nodeInfo(1,end)=10;%NE corner
nodeInfo(indexTopEnd,indexLeft)=14;%INE corner

nodeInfo(end,1)=11;%SW corner
nodeInfo(end,indexRight)=11;%SW corner
nodeInfo(indexBottomStart+1,indexRight)=15;%ISW corner

nodeInfo(end,indexLeft)=12;%SE corner
nodeInfo(end,end)=12;%SE corner
nodeInfo(indexBottomStart+1,indexLeft)=16;%ISE corner

boundOrientation=[1 1;2 2;3 3;4 4;5 2;6 3;7 1;8 4;9 5;10 6;11 7;12 8;13 9;14 10;15 11;16 12];
end