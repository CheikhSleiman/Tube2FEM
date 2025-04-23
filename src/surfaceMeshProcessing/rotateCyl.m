function hm = rotateCyl(hm,lineVec,lineVecXYperp,startNodeCoor,ThetaXInDegrees,ThetaZInDegrees) 

% Rotate cylinder around y-axis by 90
rotate(hm, [0 1 0],90,[startNodeCoor(1) startNodeCoor(2) startNodeCoor(3)])
% Other rotation 
if lineVec(1) > 0 && lineVec(2) > 0 && lineVec(3) == 0
%Rotate cylinder around z-axis by ThetaXInDegrees
rotate(hm, [0 0 1],ThetaXInDegrees,[startNodeCoor(1) startNodeCoor(2) startNodeCoor(3)])
elseif lineVec(1) < 0 && lineVec(2) > 0 && lineVec(3) == 0
%Rotate cylinder around z-axis by ThetaXInDegrees
rotate(hm, [0 0 1],ThetaXInDegrees,[startNodeCoor(1) startNodeCoor(2) startNodeCoor(3)])
elseif lineVec(1) > 0 && lineVec(2) < 0 && lineVec(3) == 0
%Rotate cylinder around z-axis by ThetaXInDegrees
rotate(hm, [0 0 1],-ThetaXInDegrees,[startNodeCoor(1) startNodeCoor(2) startNodeCoor(3)])
elseif lineVec(1) < 0 && lineVec(2) < 0 && lineVec(3) == 0
%Rotate cylinder around z-axis by ThetaXInDegrees
rotate(hm, [0 0 1],-ThetaXInDegrees,[startNodeCoor(1) startNodeCoor(2) startNodeCoor(3)])
elseif lineVec(1) == 0 && lineVec(2) > 0 && lineVec(3) == 0
%Rotate cylinder around z-axis by ThetaXInDegrees
rotate(hm, [0 0 1],ThetaXInDegrees,[startNodeCoor(1) startNodeCoor(2) startNodeCoor(3)])
elseif lineVec(1) > 0 && lineVec(2) > 0 && lineVec(3) > 0
%Rotate cylinder around z-axis by ThetaXInDegrees
rotate(hm, [0 0 1],ThetaXInDegrees,[startNodeCoor(1) startNodeCoor(2) startNodeCoor(3)])
%Rotate cylinder around lineVecXYperp by ThetaZInDegrees
rotate(hm,lineVecXYperp,ThetaZInDegrees-90,[startNodeCoor(1) startNodeCoor(2) startNodeCoor(3)])
elseif lineVec(1) < 0 && lineVec(2) > 0 && lineVec(3) > 0
%Rotate cylinder around z-axis by ThetaXInDegrees
rotate(hm, [0 0 1],ThetaXInDegrees,[startNodeCoor(1) startNodeCoor(2) startNodeCoor(3)])
%Rotate cylinder around lineVecXYperp by ThetaZInDegrees
rotate(hm,lineVecXYperp,90-ThetaZInDegrees,[startNodeCoor(1) startNodeCoor(2) startNodeCoor(3)])
elseif lineVec(1) > 0 && lineVec(2) < 0 && lineVec(3) > 0
%Rotate cylinder around z-axis by ThetaXInDegrees
rotate(hm, [0 0 1],-ThetaXInDegrees,[startNodeCoor(1) startNodeCoor(2) startNodeCoor(3)])
%Rotate cylinder around lineVecXYperp by ThetaZInDegrees
rotate(hm,lineVecXYperp,ThetaZInDegrees-90,[startNodeCoor(1) startNodeCoor(2) startNodeCoor(3)])   
elseif  lineVec(1) > 0 && lineVec(2) > 0 && lineVec(3) < 0
%Rotate cylinder around z-axis by ThetaXInDegrees
rotate(hm, [0 0 1],ThetaXInDegrees,[startNodeCoor(1) startNodeCoor(2) startNodeCoor(3)])
%Rotate cylinder around lineVecXYperp by ThetaZInDegrees
rotate(hm,lineVecXYperp,ThetaZInDegrees-90,[startNodeCoor(1) startNodeCoor(2) startNodeCoor(3)])   
elseif  lineVec(1) < 0 && lineVec(2) < 0 && lineVec(3) > 0
%Rotate cylinder around z-axis by ThetaXInDegrees
rotate(hm, [0 0 1],-ThetaXInDegrees,[startNodeCoor(1) startNodeCoor(2) startNodeCoor(3)])
%Rotate cylinder around lineVecXYperp by ThetaZInDegrees
rotate(hm,lineVecXYperp,-ThetaZInDegrees-270,[startNodeCoor(1) startNodeCoor(2) startNodeCoor(3)]) 
elseif  lineVec(1) > 0 && lineVec(2) < 0 && lineVec(3) < 0
%Rotate cylinder around z-axis by ThetaXInDegrees
rotate(hm, [0 0 1],-ThetaXInDegrees,[startNodeCoor(1) startNodeCoor(2) startNodeCoor(3)])
%Rotate cylinder around lineVecXYperp by ThetaZInDegrees
rotate(hm,lineVecXYperp,ThetaZInDegrees+270,[startNodeCoor(1) startNodeCoor(2) startNodeCoor(3)])
elseif  lineVec(1) < 0 && lineVec(2) > 0 && lineVec(3) < 0
%Rotate cylinder around z-axis by ThetaXInDegrees
rotate(hm, [0 0 1],ThetaXInDegrees,[startNodeCoor(1) startNodeCoor(2) startNodeCoor(3)])
%Rotate cylinder around lineVecXYperp by ThetaZInDegree
rotate(hm,lineVecXYperp,90-ThetaZInDegrees,[startNodeCoor(1) startNodeCoor(2) startNodeCoor(3)])
elseif  lineVec(1) < 0 && lineVec(2) < 0 && lineVec(3) < 0
%Rotate cylinder around z-axis by ThetaXInDegrees
rotate(hm, [0 0 1],-ThetaXInDegrees,[startNodeCoor(1) startNodeCoor(2) startNodeCoor(3)])
%Rotate cylinder around lineVecXYperp by ThetaZInDegrees
rotate(hm,lineVecXYperp,90-ThetaZInDegrees,[startNodeCoor(1) startNodeCoor(2) startNodeCoor(3)])
end      