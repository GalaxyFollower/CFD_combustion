function [ A,B ] = createSystemMatrices(dimX,dimY,X,Y,nodeInfo,boundOrientation , C)

% Boundary conditions 
    %0 == dirichlet with in second column the pressure (in bar); third column not important
    %1 == Neumann with in second column the pressure flux (in bar); third column not important
    %2== Robin with in second column Pinf and in third column alpha !! ONLY WEST
    %every row= boundary type: row1=west boundary         row9=Northwestcorner
    %                          row2=east boundary         row10=Northeastcorner
    %                          row3=north boundary        row11=Southwestcorner
    %                          row4=south boundary        row12=Southeastcorner
    %                          row5=east_inner boundary   %GONE row13=INW
    %                          row6=north_inner boundary  %GONE row14=INE
    %                          row7=west_inner boundary   %GONE row15=ISW
    %                          row8=south_inner boundary  %GONE row16=ISE
    
boundParameters = [1  0 0;
                    2 0 10^10;
                    1 0 0;
                    1 0 0;
                    1 0 0;
                    1 0 0;
                    1 0 0;
                    1 0 0;
                    1 0 0;
                    1 0 0;
                    1 0 0;
                    1 0 0];
                    %1 0 0;
                    %1 0 0;
                    %1 0 0;
                    %1 0 0];

%Initialising helper variables
ndof = dimY*dimX;

%Initialising system matrices
A = zeros(ndof,ndof);
B = zeros(dimX,dimY);

%Main loops for assembling the system matrices
for i = 1:dimY
    for j = 1:dimX
        
        %Helper variables for this node
        nodeType = nodeInfo(i,j);
        
        % Fill the system matrix and the right-hand side for node (i,j)
        [output1,output2]=stamp(i, j, dimX, dimY,X,Y,boundParameters,nodeType,boundOrientation,C);
        A(index(i,j,dimY), :) = output1;%row vector
        B(i,j)=output2;%B has the grid-dimension
    end
end

end

