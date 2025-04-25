function [onEdgeInput,onEdgeInputCoor,labelCount] = onEdgeInputDetection(data)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

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
ind = find(sum(onEdgeInputCoor,2)==0) ;
onEdgeInputCoor(ind,:) = [] ;




end