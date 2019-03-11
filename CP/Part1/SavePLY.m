% Function
function SavePLY(pointCloud,filename)
disp('Saving as .PLY File')
tic
    pointNum = size(pointCloud,1);
    file=fopen(filename,'w');
    fprintf(file, ['ply', '\n']);
    fprintf(file, ['format ascii 1.0\n']);
    fprintf(file, ['element vertex ', num2str(pointNum), '\n']);
    fprintf(file, ['property float x', '\n']);
    fprintf(file, ['property float y', '\n']);
    fprintf(file, ['property float z', '\n']);
    fprintf(file, ['end_header', '\n']);

    for p=1:pointNum
        % Save x, y, z coordinate as float
        fprintf(file, '%f %f %f\n', pointCloud(p,1), pointCloud(p,2), pointCloud(p,3)); 
    end
    fclose(file);
toc
end

