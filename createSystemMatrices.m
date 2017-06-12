function [ A,B ] = createSystemMatrices(dimX,dimY,X,Y,nodeInfo,boundOrientation, boundParameters )

%Initialising helper variables
ndof = dimY*dimX;
lambda = 1;
alpha = 1;

%Initialising system matrices
A = zeros(ndof,ndof);
B = zeros(dimX,dimY);

%Main loops for assembling the system matrices
for i = 1:dimY
    for j = 1:dimX
        
        %Helper variables for this node
        nodeType = nodeInfo(i,j);
        
        % Fill the system matrix and the right-hand side for node (i,j)
        [output1,output2]=stamp(i, j, dimX, dimY,X,Y,boundParameters,nodeType,boundOrientation);
        A(index(i,j,dimX,dimY), :) = output1;%row vector
        B(i,j)=output2;%B has the grid-dimension
    end
end

end

