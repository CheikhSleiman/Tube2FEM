function hm = rotateCyl(hm,lineVec,lineVecXYperp,startNodeCoor,ThetaXInDegrees,ThetaZInDegrees) 
% rotateCyl is a function used to rotate a cylinder object created in Matlab.
% 
%
%  hm = rotateCyl(hm,lineVec,lineVecXYperp,startNodeCoor,ThetaXInDegrees,ThetaZInDegrees)
%  creates a new cylinder object 'hm' after rotating a cylindrical object
%  'hm' given its centerline 'lineVec', the perpendicular to its projection
%  on the XY plane 'lineVecXYperp', and angles on projections with respect
%  to the X and Z axis 'ThetaXInDegrees' and 'ThetaZInDegrees'.
%
%   INPUTS:
%       - hm - struct, cylinder object 
%       - lineVec - vector, centerline of the cylinder
%       - lineVecXYperp - vector, perpendicular line to the projection of
%       lineVec on the XY plane
%       - startNodeCoor - vector, start node for the centerline
%       - ThetaXInDegrees - float, angle in degrees between lineVecXY and
%       [1 0 0]
%       - ThetaZInDegrees - float, angle in degrees between lineVec and
%       [0 0 1]
%   OUTPUT:
%       hm - struct, rotated cylinder object
%       
%
% -----------------------------------------------------------------------%
% Matlab seems to be inconsistant in the way it computes angles in 3D, I didn't
% find a generic way of computing angles in all the 8 quadrants and for 
% special cases so I finally computed it in the way below (& verified it case by case)



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