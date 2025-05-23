function [onEdgeInput,onEdgeInputCoor,labelCount] = onEdgeInputDetection(data)
% onEdgeInputDetection detects all input edges in a tree (skeleton or synthetic)
% represented by a data matric
%
%  [onEdgeInput,onEdgeInputCoor,labelCount] = onEdgeInputDetection(data)
%  reads 'data' which is the Tree save in 'Secomb' format and compute the
%  the input edge node index 'onEdgeInput', its coordinates and other data
% 'onEdgeInputCoor'  and gives each a label 'labelCount'.
%
%
%   INPUTS:
%       data - matrix, tree in Secomb format 
%
%   OUTPUT:
%       onEdgeInput - column vector, input edge node index
%       onEdgeInputCoor - matrix, input edge nodes coordinate
%       labelCount - column vector, input edge node label
%       
%
% -----------------------------------------------------------------------%
% Check if node index in startNode vector is not repeated in endNode vector
% Check the count for all nodes in startNode
m =size(data,1);
for i=1:m
     startNode1 = data(i,3);
     count = 0;
     for j=1:m
      endNode1   = data(j,4);
      if startNode1 == endNode1
         count = count+1;
      end
     end
     onEdgeInput(i) = count;
end

% Create onEdgeInputCoor [mx10] that contains data from Secomb data (radius, length, index)
% + coordinate of the edge nodes and that of the connected to it.
labelCount = 1;
for i =1:m
    if onEdgeInput(i) == 0 && numel(find(data(:,3)==data(i,3)))==1
    labelCount = labelCount+1;
    onEdgeInputCoor(i,1)=data(i,3);
    onEdgeInputCoor(i,2)=data(i,4);
    onEdgeInputCoor(i,3)=data(i,7); %input x
    onEdgeInputCoor(i,4)=data(i,8); %input y
    onEdgeInputCoor(i,5)=data(i,9); %input z
    onEdgeInputCoor(i,6)=data(i,10);% same branch node x
    onEdgeInputCoor(i,7)=data(i,11);% same branch node y
    onEdgeInputCoor(i,8)=data(i,12);% same branch node z
    onEdgeInputCoor(i,9)= data(i,5)/2;%radius
    onEdgeInputCoor(i,10)= labelCount;
    % scatter3(data(i,7),data(i,8), data(i,9),'MarkerEdgeColor','k',...
    %     'MarkerFaceColor',[1. 1. 1.])
    hold on
    end
end
% if count not equal zero (not edge) remove from matrix
ind = find(sum(onEdgeInputCoor,2)==0) ;
onEdgeInputCoor(ind,:) = [] ;


end