function [ Co, Neighbours ] = SetUpMesh(dimX,dimY,h,l1,l2,c2)
%Set up Mesh with X and Y points
dimY_func=dimY*ones(1,dimX);
x = 0:(l1+l2+c2)/(dimX-1):l1+l2+c2;
y1 =flip(0:h/(dimY-1):h);
size(y1)
%c1=4*h/(dimY-1) %c1=multiple of h/(dimY-1);dimY>2*multiple
%y2=c1:h/(dimY-1):h-c1;

Co=[];
Neighbours=[];

for i=1:dimX
    for j=1:dimY_func(i)
        Co=[Co; x(i),y1(j)];        
        Neighbours=[Neighbours;index(j-1,i,dimX,dimY),index(j+1,i,dimX,dimY),index(j,i-1,dimX,dimY),index(j,i+1,dimX,dimY)];
    end
end
end

