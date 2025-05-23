function [GmshStatus] = writeGmsh(fileName,meshOutput)
% writeGmsh is a function that write .msh format files from TetGen output
% 
%
% [GmshStatus] = writeGmsh(fileName,meshOutput) takes a gmsh file name 
% 'fileName'as input and a tetGen mesh struct and produces a version 2
% ASCII Gmsh file format.
%
%   INPUTS:
%       - fileName - string, Gmsh file name that you want to write 
%       - meshOutput - struct, TetGen mesh
%   OUTPUT:
%       Gmsh file format (.msh)
%       
%
% -----------------------------------------------------------------------%


% Read mesh boundaries (surface mesh and interfaces if any) and labels
triangleMat = [meshOutput.facesBoundary meshOutput.boundaryMarker];
triangleMatSort = sortrows(triangleMat,4);

% Read tet elements and elementMaterialID
tetMat = [meshOutput.elements meshOutput.elementMaterialID];
tetMat = sortrows(tetMat,5);

% open file in writing mode
fid = fopen(fileName,'wt');
% Specify this line to say it is a version 2 ASCII format
fprintf(fid, '$MeshFormat \n2.2 0 8 \n$EndMeshFormat\n') ;
% Physical names are labels for surfaces and volumes
fprintf(fid,'\n$PhysicalNames\n');
% Labelled volumes
labelledVols = max(meshOutput.elementMaterialID);
% Labelled surfaces
labelledSurfs = max(meshOutput.boundaryMarker);
% Their total number
fprintf(fid, '%d \n',labelledVols+labelledSurfs);

% Surface physical number
for i=1:max(meshOutput.boundaryMarker)
    fprintf(fid,'%d %d "%d"\n',2,i,i);
end

% Volume physical number
for i=1:max(meshOutput.elementMaterialID)
    fprintf(fid,'%d %d "%d"\n',3,i,i);
end
fprintf(fid, '$EndPhysicalNames \n');


% Nodes 
fprintf(fid, '$Nodes\n');
fprintf(fid,'%d\n',size(meshOutput.nodes,1));

for i =1:size(meshOutput.nodes,1)
      fprintf(fid, '%d %.16f %.16f %.16f\n',i,meshOutput.nodes(i,1),meshOutput.nodes(i,2),meshOutput.nodes(i,3));
end

% Elements
fprintf(fid, '$EndNodes \n$Elements\n');
fprintf(fid,'%d\n',size(meshOutput.facesBoundary,1)+size(meshOutput.elements,1));

%Triangles
for i =1:size(meshOutput.boundaryMarker,1)
    fprintf(fid, '%d %d %d %d %d %d %d %d\n',i,2,2,triangleMatSort(i,4),triangleMatSort(i,4),triangleMatSort(i,1),triangleMatSort(i,2),triangleMatSort(i,3));
end

%Tets
for i =1:size(meshOutput.elements,1)
    fprintf(fid, '%d %d %d %d %d %d %d %d %d\n',size(meshOutput.boundaryMarker,1)+i,4,2,tetMat(i,5),tetMat(i,5),tetMat(i,1),tetMat(i,2),tetMat(i,3),tetMat(i,4));
end

fprintf(fid, '$EndElements\n');
fclose(fid);

GmshStatus = fprintf('Gmsh Version 2 ASCII written\n');
end

