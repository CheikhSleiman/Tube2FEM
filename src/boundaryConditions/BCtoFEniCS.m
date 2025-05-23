function BCtoFEniCS(meshOutput,Fc,Vc)
% BCtoFEniCS Generates a 'BC.txt' where all BCs are written in a way that is
% compatible with FEniCS
%
%   BCtoFEniCS(meshOutput, Fc, Vc) reads the TetGen volumetric 
%   mesh struct 'meshOutput' and surface mesh nodes 'Vc' and faces 'FC' 
%   and writes a text file containing a FEniCS compatible BC format.
%   
%   INPUTS:
%       meshOutput  - struct, meshOutput from TetGen
%       Fc - matrix, surface mesh faces
%       Vc - matrix, surface mesh nodes
%
%   OUTPUT:
%       Textfile containing FEniCS compatible BC syntax
%
% -----------------------------------------------------------------------%

% Get group counts based on the labels of the surface mesh
[GC,GR] = groupcounts(meshOutput.boundaryMarker) ;

% Get the two largest groups - the largest labelled surface is the wall
% and the second largest is the inlet surface. This only applies for
% hierarchical trees
twoLargest = maxk(GC,2);

% Wall surface
firstMax  = twoLargest(1);
% Inlet surface 
secondMax = twoLargest(2);
% Get indices
indWall= find(GC== firstMax);
indu = find(GC == secondMax);

%find normal to surf indu (inlet)
indSurfindu = find(meshOutput.boundaryMarker == indu);
Nu=abs(round(mean(patchNormal(Fc(indSurfindu,1:3),Vc))));
ind1 = find(Nu==1);

labelledSurfs = max(meshOutput.boundaryMarker);

% Open text file
fid1 = fopen('BC.txt','wt');
% Write the velocity expression applied on the inlet
fprintf(fid1,'u_in = Expression(''-0.025*t'', degree = 1, t=0)\n\n');
% Write the outlet pressure expression
fprintf(fid1,'p_out = Expression(''0*t'', degree = 1, t=0)\n\n');
% Loop over all outlets and assign p_out for all its labels
for i =1:max(labelledSurfs)
    if i ~= indu && i ~= indWall
        fprintf(fid1, 'bcp_outflow%d = DirichletBC(Q, p_out, surf_markers, %d)\n',i,i); 
    end
end

fprintf(fid1, '\n');
% Assign null velocity for walls
fprintf(fid1, 'bcu_walls     = DirichletBC(V, Constant((0, 0, 0)),surf_markers, %d)\n',indWall);
% Assign u_in at the normal direction to the inlet surface
fprintf(fid1, 'bcu_inflow1    = DirichletBC(V.sub(%d), u_in , surf_markers, %d)',ind1-1,indu);
fprintf(fid1, '\n');
fprintf(fid1, '%s','bcu = [bcu_walls,bcu_inflow1]');
fprintf(fid1, '\n');

% Write the velocity bcu and pressure bcp vector that assembles everything
indAll = 1:max(labelledSurfs);
induAll = [indu, indWall];
flag = ~ismember(indAll,induAll);
index = find(flag);
maxp = max(index);
fprintf(fid1, 'bcp =[');
for i =1:max(labelledSurfs)
    if i ~= indu && i ~= indWall && i~= maxp
        fprintf(fid1, 'bcp_outflow%d,',i);
    elseif i ~= indu && i ~= indWall && i== maxp
        fprintf(fid1, 'bcp_outflow%d',i);
    end
end
fprintf(fid1, ']');
% Close the text file!
fclose(fid1);

fprintf('BCs are written correctly in BC.txt \n');
end

