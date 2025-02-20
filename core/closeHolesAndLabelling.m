function [Fc,Vc,Cc] = closeHolesAndLabelling(Fc,Vc,Cc)
Ec=patchBoundary(Fc,Vc);

optionStruct.outputType='label';
G=tesgroup(Ec,optionStruct);
i=0
while 1
    
    Ec=patchBoundary(Fc,Vc);

    optionStruct.outputType='label';
    G=tesgroup(Ec,optionStruct);
    Ec=Ec(G==1,:);
    indCut=edgeListToCurve(Ec);

    %---------- New

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


    % -------
    
    % regiontrimesh3D
    regionCell={Vc(indCut(1,1:end-1),:)};
    
    pointSpacing=mean(patchEdgeLengths(Fc,Vc)); %Desired point spacing
    resampleCurveOpt=0;
    interpMethod='linear'; %or 'natural'
    [Fsi,Vsi]=regionTriMesh3D(regionCell,pointSpacing,resampleCurveOpt,interpMethod);
    gpatch(Fsi,Vsi,'y','k',1)
    axisGeom
    hold on
    Csi = ones(size(Fsi,1),1);
    % join element set
    [Fc,Vc,Cc]=joinElementSets({Fc,Fsi},{Vc,Vsi},{Cc,Csi+max(Cc(:))});
    [Fc,Vc]=mergeVertices(Fc,Vc);
    i=i+1
    if max(G)==1
        break
    end
end
% labelled surfaces should at least have 10 elements, orelse they are
% defects and therefore joined to wall surface
[GC,GR] = groupcounts(Cc);
for i=1:size(GC,1)
    if GC(i) <10
        indGC = find(Cc==i);
        Cc(indGC)=ones(size(indGC,1),1);
    end
end

[GC,GR] = groupcounts(Cc);
for i = 1:size(GR,1)
    indGR = find(Cc==GR(i));
    Cc(indGR)=i*ones(size(indGR,1),1);
end

% patch('faces',F,'vertices',V,'FaceColor','flat','CData',C,'FaceAlpha',0.2,'EdgeColor','none');

end