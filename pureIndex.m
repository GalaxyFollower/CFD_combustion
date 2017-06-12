function [ index ] = pureIndex( i,j,dimY )
%PUREINDEX Summary of this function goes here
%   Detailed explanation goes here

    index = i + (j-1) .* dimY;
end

