function ExportBCs(FEFolder,Vsurf,oxygenSurf)



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

% Close the file
fclose(fid1);
fclose(fid2);