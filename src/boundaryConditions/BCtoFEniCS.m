function BCtoFEniCS(meshOutput,Fc,Vc)

[GC,GR] = groupcounts(meshOutput.boundaryMarker) ;
twoLargest = maxk(GC,2);

firstMax  = twoLargest(1);
secondMax = twoLargest(2);
indWall= find(GC== firstMax);
indu = find(GC == secondMax);

%find normal to surf indu
indSurfindu = find(meshOutput.boundaryMarker == indu);
Nu=abs(round(mean(patchNormal(Fc(indSurfindu,1:3),Vc))));
ind1 = find(Nu==1);

labelledSurfs = max(meshOutput.boundaryMarker);

fid1 = fopen('BC.txt','wt');
fprintf(fid1,'u_in = Expression(''-0.025*t'', degree = 1, t=0)\n\n');
fprintf(fid1,'p_out = Expression(''0*t'', degree = 1, t=0)\n\n');
for i =1:max(labelledSurfs)
    if i ~= indu && i ~= indWall
        fprintf(fid1, 'bcp_outflow%d = DirichletBC(Q, p_out, surf_markers, %d)\n',i,i); 
    end
end

fprintf(fid1, '\n');
fprintf(fid1, 'bcu_walls     = DirichletBC(V, Constant((0, 0, 0)),surf_markers, %d)\n',indWall);
fprintf(fid1, 'bcu_inflow1    = DirichletBC(V.sub(%d), u_in , surf_markers, %d)',ind1-1,indu);
fprintf(fid1, '\n');
fprintf(fid1, '%s','bcu = [bcu_walls,bcu_inflow1]');
fprintf(fid1, '\n');

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
fclose(fid1);

fprintf('BCs are written correctly in BC.txt \n');
end

