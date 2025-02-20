close all force
clear all
clc

%% Adding path to Tube2FEM core directory
CurrentFolder = pwd;
TopFolder = fileparts(pwd);
TopTopFolder = fileparts(fileparts(pwd));
coreFolder = strcat(TopTopFolder,'\core');
addpath(coreFolder);

%% Read Tif
VolTif = tiffreadVolume('Input\RootSoil_tomo.tif');
V = volshow(fliplr(VolTif));
V.Parent.BackgroundColor = [1 1 1];
V.Parent.BackgroundGradient='off';
V.RenderingStyle = 'GradientOpacity';
V.Parent.LightColor = [1 1 0];
V.Parent.LightPositionMode = 'left';
V.Parent.OrientationAxes = 'off';


%% Get .stl from binary image
logicVoxels=(VolTif>0);
[Frs,Vrs]=im2patch(VolTif,logicVoxels);
[Frs,Vrs]=patch2tri(Frs,Vrs);


cFigure; 
gpatch(Frs,Vrs,'w','k');
axisGeom;
axis off
view([138,14])


%% Remesh
optionStruct1.nb_pts= 100000;
[Frr,Vrr]=ggremesh(Frs,Vrs,optionStruct1);
Vrr = Vrr*0.0002;
cFigure;
gpatch(Frr,Vrr,'y','k')
axisGeom;
axis off

%% Smooth
nSteps = 2;
cPar2.Method='HC'; %Smooth method
cPar2.n=nSteps; %Number of iterations
[Vrr]=patchSmooth(Frr,Vrr,[],cPar2);
cFigure;
gpatch(Frr,Vrr,'y','k')
axisGeom;
axis off

%% Save surface mesh as .stl
fileName= 'Mesh\smoothMesh.stl';
patch2STL(fileName,Vrr,Frr,[],'SmoothMesh');

%% Repair Mesh 
MeshPath = strcat(CurrentFolder,'\Mesh');
MeshPath = strrep(MeshPath, '\', '/');

pyrunfile("RepairSurfaceMesh.py",Path=MeshPath)
cd Mesh
fileName= 'repairedMesh.stl';
[stlStruct] = import_STL(fileName);
F=stlStruct.solidFaces{1}; %Faces
V=stlStruct.solidVertices{1}; %Vertices
[Frr,Vrr]=mergeVertices(F,V); % Merging nodes

%% Check if watertight after repairing
mesh = surfaceMesh(Vrr,Frr);
TF = isWatertight(mesh);

%% Slicing
Pcut = [0.04,0.04,0.15];
C=[]; 
snapTolerance=mean(patchEdgeLengths(F,V))/100;
n = [0,0,1]
[Frr,Vrr,~,logicSide]=triSurfSlice(Frr,Vrr,C,Pcut,n,snapTolerance);
Frr = Frr(logicSide,:);
cFigure;
gpatch(Frr,Vrr,'y','k')
axisGeom;
axis off

%% Remesh open surface mesh
nb_pts = size(Vrr,1);
optionStruct2.nb_pts = nb_pts;
[Frr,Vrr] = ggremesh(Frr,Vrr,optionStruct2);
cFigure;
gpatch(Frr,Vrr,'w','k');
axisGeom;
axis off

%% Close holes and labelling
cFigure
Crr = ones(size(Frr,1),1);
[Frr,Vrr,Crr]=closeHolesAndLabelling(Frr,Vrr,Crr);
Crr(Crr==1)=-1;
Crr(Crr==2)=1;
Crr(Crr==-1)=2;
gpatch(Frr,Vrr,Crr,'k');
axisGeom
axis off

%% Create Box
boxDim=[0.06 0.08 0.15]; %Width in each direction
pointSpacing=0.005; %Desired point spacing
[Fbox,Vbox,Cbox]=triBox(boxDim,pointSpacing);

% Move Box
Vx = 0.05;
Vy = 0.04;
Vz = 0.078;
MoveBox = ones(size(Vbox,1),size(Vbox,2));
MoveBox = [Vx*MoveBox(:,1),Vy*MoveBox(:,2),Vz*MoveBox(:,3)];
Vbox = Vbox+MoveBox;

% Visualisation
cFigure; 
% Crr = ones(size(Frr,1),1);
patch('faces',Frr,'vertices',Vrr,'FaceColor','flat','CData',Crr,'FaceAlpha',0.8,'EdgeColor','none');
axisGeom;
% Cbox = 2*ones(size(Fbox,1),1);
patch('faces',Fbox,'vertices',Vbox,'FaceColor','flat','CData',Cbox,'FaceAlpha',0.2,'EdgeColor','w');
axis off


%% Join surface sets
Vnet = Vrr;
Fnet = Frr;
Cnet = Crr;
[F,V,C]=joinElementSets({Fbox,Fnet},{Vbox,Vnet},{Cbox,Cnet});
C(1:size(Cbox(:),1))=Cbox;
C(size(Cbox(:),1)+1:size(C,1)) = 6*ones((size(C,1)-size(Cbox,1)),1)+Cnet;

% Find interior points
% [V_region1]=getInnerPoint({Fbox,Fnet},{Vbox,Vnet});
[V_region1] = [0.03 0.01 0.03];
% [V_region2]=getInnerPoint(Fnet,Vnet);
[V_region2] = [0.04 0.04 0.15];
% [V_region2] = [V_region2; [475e-6 288e-6 180e-6]];

V_regions= [V_region1;V_region2];

% Volume parameters
[vol1]=tetVolMeanEst(Fbox,Vbox);
[vol2]=tetVolMeanEst(Fnet,Vnet);

regionTetVolumes=[vol1;vol2]; %Element volume settings
stringOpt='-pq1.2AaY'; %Tetgen options
modelName = 'test';

% Mesh inputs
% Create tetgen input structure
inputStruct.stringOpt=stringOpt; %Tetgen options
inputStruct.Faces=F; %Boundary faces
inputStruct.Nodes=V; %Nodes of boundary
inputStruct.faceBoundaryMarker=C;
inputStruct.regionPoints=V_regions; %Interior points for regions
inputStruct.regionA=regionTetVolumes; %Desired tetrahedral volume for each region
inputStruct.modelName = modelName;

% Mesh model using tetrahedral elements using tetGen
[meshOutput]=runTetGen(inputStruct); %Run tetGen

% Mesh Output
E=meshOutput.elements; %The elements
V=meshOutput.nodes; %The vertices or nodes
CE=meshOutput.elementMaterialID; %Element material or region id
Fb=meshOutput.facesBoundary; %The boundary faces
Cb=meshOutput.boundaryMarker; %The boundary markers

for i=1:size(meshOutput.elementMaterialID,1)
    if meshOutput.elementMaterialID(i)==1
        meshOutput.elementMaterialID(i) = 2;
    elseif meshOutput.elementMaterialID(i)==-2
        meshOutput.elementMaterialID(i) = 1;
    end
end

% Visualization
hf=cFigure; hold on;

% Visualizing using |meshView|
optionStruct.hFig=hf;

meshView(meshOutput,optionStruct);

colormap(matplotlibColormap(3,2))
axis off
view([-37,21])
camlight('left')
% print(gcf,'RootSoilEmbeddedMesh.png','-dpng','-r900');
% print(gcf,'RootSoilEmbeddedMeshZoom2.png','-dpng','-r1500');

%% Write version 2 ASCII .msh mesh format 
GmshFileName = 'RootSoilGmshFormat.msh';
[status] = writeGmsh(GmshFileName,meshOutput);

%% Convert .msh to .xdmf via meshio
meshPath = strcat(CurrentFolder,'\Mesh');
meshPath = strrep(meshPath, '\', '/');
% Split the path into parts
pathParts = strsplit(meshPath, '/');
% Reconstruct the path without the first directory
newMeshPath = strjoin(pathParts(2:end), '/'); 
rootPath = "../../../../../../../../../../../../../../../../../mnt/c/";
wslMeshPath = strcat(rootPath,newMeshPath,"/meshioConvertMesh.py");
wslMeshPath= sprintf('"%s"', wslMeshPath);

fid1 = fopen('meshInfo.txt','wt');
meshPath1 = sprintf('"%s"', strcat(rootPath,newMeshPath));
fprintf(fid1,strcat('meshPath=',meshPath1));
fprintf(fid1,'\n');
GmshFileName= sprintf('"%s"', GmshFileName);
fprintf(fid1,strcat('GmshFileName=',GmshFileName));
fclose(fid1);

script1 = fileread('meshInfo.txt');
script2 = fileread('../../../core/ConvertGmshToXdmf.py');

fid2 = fopen('meshioConvertMesh.py','wt');
fprintf(fid2,script1);
fprintf(fid2,'\n');
fprintf(fid2,script2);
fprintf(fid2,'\n');
fclose(fid2);

wslPath = '"C:\Windows\System32\wsl.exe"';
runString = strcat(wslPath,' python3',append(' ',wslMeshPath));
[runStatus,runOut]=system(runString,'-echo');

