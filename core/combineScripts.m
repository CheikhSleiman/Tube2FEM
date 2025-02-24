function combineScripts
scriptBeforeBC = fileread('NavierStokesTimeDependent_beforeBC.py');
bcContent = fileread('BC.txt');
scriptAfterBC = fileread('NavierStokesTimeDependent_afterBC.py');


fid1 = fopen('combinedNavierStokes.py','wt');
fprintf(fid1,scriptBeforeBC);
fprintf(fid1,'\n');
fprintf(fid1,bcContent);
fprintf(fid1,'\n');
fprintf(fid1, '%s',scriptAfterBC);
fclose(fid1);
end