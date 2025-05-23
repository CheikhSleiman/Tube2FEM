function [Fc,Vc,Cc] = closeHolesAndLabelling(Fc,Vc,Cc)
% closeHolesAndLabelling closes all the open holes done by the slicing 
% operation by creating new surfaces and assigning unique labels to each.
%
%   [Fc,Vc,Cc] = closeHolesAndLabelling(Fc,Vc,Cc) reads the open surface mesh
%  faces 'Fc', vertices 'Vc' and labels 'Cc' (initially all labels are equal 1)
%  and create a new water tight labelled surface mesh (Fc,Vc,Cc)
%
%
%   INPUTS:
%       Fc - matrix, surface mesh faces
%       Vc - matrix, surface mesh nodes
%       Cc - column vector - surface mesh labels
%
%   OUTPUT:
%       Fc - matrix, surface mesh faces (with new faces for inlets/outlets)
%       Vc - matrix, surface mesh nodes 
%       Cc - column vector - surface mesh labels (with labels for inlets/outlets)
%       
%
% -----------------------------------------------------------------------%


% Get the boundary of the open surface mesh
Ec=patchBoundary(Fc,Vc);
% Sort boundary elements by label
optionStruct.outputType='label';
G=tesgroup(Ec,optionStruct);
i=0;
% As long as there is open surface in the mesh continue looping
while 1
    
    Ec=patchBoundary(Fc,Vc);

    optionStruct.outputType='label';
    G=tesgroup(Ec,optionStruct);
    % Get the first labelled group (Ec will be decreasing in size after each iteration)
    Ec=Ec(G==1,:);

    indCut=edgeListToCurve(Ec);

    % sometimes the nodes on the boundary are not coplanar, this piece of
    % code will fix this issue:

    Vreg = Vc(indCut(1,1:end-1),:);
    Creg = ones(size(Vreg,1),1);
    meanV = mean(Vreg);
    %center data
    Vregc = Vreg - meanV;
    % Perform PCA
    [coeff, ~, latent] = pca(Vregc);

    % The normal vector is the principal component corresponding to the smallest eigenvalue
    [~, min_idx] = min(latent);
    normal_vector = coeff(:, min_idx);
    normal_vector = normal_vector';
    normal_vector = normal_vector / norm(normal_vector);
    Vreg_p = zeros(size(Vreg));
    for i = 1:size(Vreg, 1)
        Vregi = Vreg(i, :);
        vector = Vregi - meanV;
        distance = dot(vector, normal_vector);
        Vreg_p(i, :) = Vregi - distance * normal_vector;
    end
    Vc(indCut(1,1:end-1),:) = Vreg_p;


    % Now that the nodes are coplanar, 3D mesh the surface
    
    regionCell={Vc(indCut(1,1:end-1),:)};
    % regionTriMesh3D
    pointSpacing=mean(patchEdgeLengths(Fc,Vc)); %Desired point spacing
    resampleCurveOpt=0;
    interpMethod='linear'; %or 'natural'
    [Fsi,Vsi]=regionTriMesh3D(regionCell,pointSpacing,resampleCurveOpt,interpMethod);
    % plot
    gpatch(Fsi,Vsi,'y','k',1)
    axisGeom
    hold on
    % use same labels for the plot
    Csi = ones(size(Fsi,1),1);
    
    % join element set (the newly formed meshed inlet/outlet surface + mesh wall)
    % Add a distinct label to each inlet/outlet
    [Fc,Vc,Cc]=joinElementSets({Fc,Fsi},{Vc,Vsi},{Cc,Csi+max(Cc(:))});
    [Fc,Vc]=mergeVertices(Fc,Vc);
    i=i+1;
    % Break the while loop when we have one group remaining
    if max(G)==1
        break
    end
end
end