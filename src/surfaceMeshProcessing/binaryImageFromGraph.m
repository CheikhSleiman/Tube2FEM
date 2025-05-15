function [ImBinary] = binaryImageFromGraph(data,voxelSize)
% binaryImageFromGraph generate a binary image from a synthetic or
% skeletonised graph (presented by the secomb format)
%
%  [ImBinary] = binaryImageFromGraph(data) reads 'data' which is the graph
%  saved in 'Secomb' format and generate 'Imbinary' a 3D binary image. This
%  could be seen as an inverse skeletonisation function
%
%   INPUTS:
%       data - matrix, graph in Secomb format 
%       voxelSize - float, voxel size value
%   OUTPUT:
%       ImBinary - 3D matrix, 3D binary image of the graph
%       
%
% -----------------------------------------------------------------------%



a = linspace(0,1);
b = zeros(1,100);
c = zeros(1,100);
myMap = [a' b' c'];

maxX = 10*(ceil(max(data(:,7))/10.)+1);
maxY = 10*(ceil(max(data(:,8))/10.)+1);
maxZ = 10*(ceil(max(data(:,9))/10.)+1);


ImBinary = zeros(maxX,maxY,maxZ);

x = [0:voxelSize:size(ImBinary,1)];
y = [0:voxelSize:size(ImBinary,2)];
z = [0:voxelSize:size(ImBinary,3)];

cFigure

for i = 1:size(data,1)
     startNode = data(i,3);
     endNode  = data(i,4);
     barRadius = data(i,5);
     barLength = data(i,6);
     startNodeCoor = [data(i,7) data(i,8) data(i,9)];
     endNodeCoor = [data(i,10) data(i,11) data(i,12)];
     line([startNodeCoor(1) endNodeCoor(1)],[startNodeCoor(2) endNodeCoor(2)],[startNodeCoor(3) endNodeCoor(3)])     


     lineVec = endNodeCoor - startNodeCoor;
     quiver3(startNodeCoor(1),startNodeCoor(2),startNodeCoor(3),lineVec(1),lineVec(2),lineVec(3),'k','LineWidth',(20*barRadius/100),'ShowArrowHead','off')
     hold on
     grid on
     lineVecXY = [lineVec(1) lineVec(2) 0];
     lineVecXZ = [lineVec(1) 0 lineVec(3)];

     %Plot perpendicular vector to lineVecXY
     if sum(lineVecXY ~= 0) == 0
        lineVecXYperp = [1 0 0];
        % ThetaXInDegrees = 0
        % ThetaZInDegrees = -90
     else
        lineVecXYperp = find_perp(lineVecXY);
        lineVecXYperp = [lineVecXYperp(1) lineVecXYperp(2) 0];
        
     end
        %Angles calculation
        vx = [1 0 0];
        CosTheta = max(min(dot(lineVecXY,vx)/(norm(lineVecXY)*norm(vx)),1),-1);
        ThetaXInDegrees = real(acosd(CosTheta));
        vz = [0 0 1];
        CosTheta = max(min(dot(lineVec,vz)/(norm(lineVec)*norm(vz)),1),-1);
        ThetaZInDegrees = real(acosd(CosTheta));
     % Create Cylinder
     [X,Y,Z] = cylinder(barRadius); % use cylinder(barRadius,12) to fasten it up
     Z = Z*barLength;
     maxRadius = max(data(:,5));
     hm = surface(X+startNodeCoor(1),Y+startNodeCoor(2),Z+startNodeCoor(3),'facecolor',barRadius/maxRadius*[1 0 0],'LineStyle','none','facealpha',0.5);
     %Rotate cylinder function
     hm = rotateCyl(hm,lineVec,lineVecXYperp,startNodeCoor,ThetaXInDegrees,ThetaZInDegrees);      

% Define a cylinder domain
xmin = min(min(hm.XData));
ymin = min(min(hm.YData));
zmin = min(min(hm.ZData));

xmax = max(max(hm.XData));
ymax = max(max(hm.YData));
zmax = max(max(hm.ZData));
      
imin = find(x < round(xmin,2)-voxelSize);
imin = max(imin);
imax = find(x > round(xmax,2)+voxelSize);
imax = min(imax);

jmin = find(y < round(ymin,2)-voxelSize);
jmin = max(jmin);
jmax = find(y > round(ymax,2)+voxelSize);
jmax = min(jmax);

kmin = find(z < round(zmin,2)-voxelSize);
kmin = max(kmin);
kmax = find(z > round(zmax,2)+voxelSize);
kmax = min(kmax);
      
% % point in cylinder: dist(P,cylinderAxis)<R + dist(Q,centerOfCylAxis)<barLength/2, Q being the projection of P on CylAxis --- Validation in PointInCyliner.m file
m1 = cross(startNodeCoor,endNodeCoor);
midPt = (startNodeCoor+endNodeCoor)/2;

threshRadius = 1;

for i =imin:imax
    for j =jmin:jmax
        for k = kmin:kmax
            pt = [x(i) y(j) z(k)];
            dist = point_to_line(pt,startNodeCoor,endNodeCoor);
            ptProj = pt + (cross(lineVec,(m1+cross(lineVec,pt))))/(norm(lineVec))^2;
            distProj = norm(ptProj-midPt);
            
            if dist < barRadius && distProj < barLength/2  && barRadius < threshRadius
                ImBinary(i,j,k) = 255;
            elseif dist < barRadius && distProj < barLength/2 && barRadius > threshRadius 
                ImBinary(i,j,k) = 255;               
            end    
        end
    end
end     


end

set(gca,'xdir','reverse','zdir','reverse')
xlim([0 maxX])
ylim([0 maxY])
zlim([0 maxZ])
xlabel('x [\mum]')
ylabel('y [\mum]')
zlabel('z [\mum]')
colormap(myMap)
minRadius = min(data(:,5));
% caxis([minRadius, maxRadius])
c = colorbar;
c.Label.String = 'Vessel diameter [\mum]';
hold on



% Spheres at the nodes
for i = 1:size(data,1)
     startNode = data(i,3);
     endNode  = data(i,4);
     barRadius = data(i,5);
     startNodeCoor = [data(i,7) data(i,8) data(i,9)];
     endNodeCoor = [data(i,10) data(i,11) data(i,12)];

      % Create Sphere
      [X,Y,Z] = sphere; %use sphere(12) to fasten it up
      r = barRadius;
      X2 = X * r;
      Y2 = Y * r;
      Z2 = Z * r;
      
      hm = surface(X2+startNodeCoor(1),Y2+startNodeCoor(2),Z2+startNodeCoor(3),'facealpha',0.5);
      hold on
      hm2 = surface(X2+endNodeCoor(1),Y2+endNodeCoor(2),Z2+endNodeCoor(3),'facealpha',0.5);


% Define a cylinder domain
      xmin = min(min(hm.XData));
      ymin = min(min(hm.YData));
      zmin = min(min(hm.ZData));

      xmax = max(max(hm.XData));
      ymax = max(max(hm.YData));
      zmax = max(max(hm.ZData));
      
      imin = find(x < round(xmin,2)-voxelSize);
      imin = max(imin);
      imax = find(x > round(xmax,2)+voxelSize);
      imax = min(imax);

      jmin = find(y < round(ymin,2)-voxelSize);
      jmin = max(jmin);
      jmax = find(y > round(ymax,2)+voxelSize);
      jmax = min(jmax);

      kmin = find(z < round(zmin,2)-voxelSize);
      kmin = max(kmin);
      kmax = find(z > round(zmax,2)+voxelSize);
      kmax = min(kmax);
      
        for i =imin:imax
            for j =jmin:jmax
                for k = kmin:kmax
                    pt = [x(i) y(j) z(k)];
                    center = [startNodeCoor(1) startNodeCoor(2) startNodeCoor(3)]; 
                    dist = norm(pt-center);
           
                    if dist < barRadius
                        ImBinary(i,j,k) = 255;               
                    end    
                end
            end
        end     


      % Define a sphere domain
      xmin2 = min(min(hm2.XData));
      ymin2 = min(min(hm2.YData));
      zmin2 = min(min(hm2.ZData));

      xmax2 = max(max(hm2.XData));
      ymax2 = max(max(hm2.YData));
      zmax2 = max(max(hm2.ZData));
      
      imin2 = find(x < round(xmin2,2)-voxelSize);
      imin2 = max(imin2);
      imax2 = find(x > round(xmax2,2)+voxelSize);
      imax2 = min(imax2);

      jmin2 = find(y < round(ymin2,2)-voxelSize);
      jmin2 = max(jmin2);
      jmax2 = find(y > round(ymax2,2)+voxelSize);
      jmax2 = min(jmax2);

      kmin2 = find(z < round(zmin2,2)-voxelSize);
      kmin2 = max(kmin2);
      kmax2 = find(z > round(zmax2,2)+voxelSize);
      kmax2 = min(kmax2);

        for i =imin2:imax2
            for j =jmin2:jmax2
                for k = kmin2:kmax2
                    pt = [x(i) y(j) z(k)];
                    center = [endNodeCoor(1) endNodeCoor(2) endNodeCoor(3)]; 
                    dist = norm(pt-center);
           
                    if dist < barRadius
                        ImBinary(i,j,k) = 255;               
                    end    
                end
            end
        end 

end



end