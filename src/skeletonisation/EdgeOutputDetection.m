function [onEdgeOutput,onEdgeOutputCoor] = onEdgeOutputDetection(data,labelCount)
% onEdgeOutputDetection detects all input edges in a tree (skeleton or synthetic)
% represented by a data matric
%
%  [onEdgeOutput,onEdgeOutputCoor,labelCount] = onEdgeOutputDetection(data)
%  reads 'data' which is the Tree save in 'Secomb' format and compute the
%  the input edge node index 'onEdgeOutput', its coordinates and other data
% 'onEdgeOutputCoor'  and gives each a label 'labelCount'.
%
%
%   INPUTS:
%       data - matrix, tree in Secomb format 
%
%   OUTPUT:
%       onEdgeOutput - column vector, input edge node index
%       onEdgeOutputCoor - matrix, input edge nodes coordinate
%       labelCount - column vector, input edge node label
%       
%
% -----------------------------------------------------------------------%
% Check if node index in endNode vector is not repeated in endNode vector
% Check the count for all nodes in endNode

m = size(data,1);
for i=1:m
     startNode1 = data(i,4);
     count = 0;
     for j=1:m
      endNode1   = data(j,3);
      if startNode1 == endNode1
         count = count+1;
      end
     end
     onEdgeOutput(i) = count;
end

% Create onEdgeOutputCoor [mx10] that contains data from Secomb data (radius, length, index)
% + coordinate of the edge nodes and that of the connected to it.
for i =1:m
    if onEdgeOutput(i) == 0 && numel(find(data(:,4)==data(i,4)))==1
    labelCount = labelCount+1;
    onEdgeOutputCoor(i,1)=data(i,4);
    onEdgeOutputCoor(i,2)=data(i,3);
    onEdgeOutputCoor(i,3)=data(i,10);   %output x
    onEdgeOutputCoor(i,4)=data(i,11);   %output y
    onEdgeOutputCoor(i,5)=data(i,12);   %output z
    onEdgeOutputCoor(i,6)=data(i,7);  %same branch x
    onEdgeOutputCoor(i,7)=data(i,8);  %same branch y
    onEdgeOutputCoor(i,8)=data(i,9);  %same branch z
    onEdgeOutputCoor(i,9)= data(i,5)/2;%radius
    onEdgeOutputCoor(i,10)=labelCount;
    % scatter3(data(i,10),data(i,11), data(i,12),'MarkerEdgeColor','k',...
    %     'MarkerFaceColor',[1. 1. 1.])
    hold on
    end
end
% if count not equal zero (not edge) remove from matrix
ind = find(sum(onEdgeOutputCoor,2)==0) ;
onEdgeOutputCoor(ind,:) = [] ;
end