function ExportBCs(FEFolder,Vsurf,oxygenSurf)
% ExportBcs Generates two text files 'BC_points.txt' and 'BC_values.txt'
% and is used for 3D-1D problems. Each nodes on the outer surface of the
% tubular object will be assigned a pointwise BC value.
%
%   ExportBCs(FEFolder,Vsurf,oxygenSurf) reads the FE directory and writes
%   in it the two text files mentioned above. 'Vsurf' is the surface mesh nodes
%   and 'oxygenSurf' is the value of the field project on that surface 
%   (in the paper it is an oxygen field).
%   
%   INPUTS:
%       FEFolder - string, FE directory
%       Vsurf - matrix, surface mesh nodes
%       oxygenSurf - column vector, surface mesh field values
%
%   OUTPUT:
%       Two textfiles containing the surface mesh nodes and surface mesh
%       nodal field values, respectively.
%
% -----------------------------------------------------------------------%


% Create the two text files
fileName1 = strcat(FEFolder,'\BC_points.txt');
fileName2 = strcat(FEFolder,'\BC_values.txt');
fid1 = fopen(fileName1,'wt');
fid2 = fopen(fileName2,'wt');

[row, col] = size(Vsurf);
for i = 1:row
    fprintf(fid2, '%d', oxygenSurf(i));  % Assuming integers in the matrix

    for j = 1:col
        if j==1 || j==2
            fprintf(fid1, '%d ', Vsurf(i, j));  % Assuming integers in the matrix
        else
            fprintf(fid1, '%d', Vsurf(i, j));  % Assuming integers in the matrix
        end
    end
    if i~=row
        fprintf(fid1, '\n');  % Newline after each row
        fprintf(fid2, '\n');  % Newline after each row
    end
end

% Close the files
fclose(fid1);
fclose(fid2);