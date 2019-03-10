% Function
function SavePLY(depthMap,filename)
disp('Saving as .PLY File')
tic
fid=fopen(filename,'w');
fprintf(fid, ['ply', '\n']);
fprintf(fid, ['format ascii 1.0\n']);
fprintf(fid, ['element vertex ', num2str(size(find(depthMap(:,:,3)),1)), '\n']);
fprintf(fid, ['property float x', '\n']);
fprintf(fid, ['property float y', '\n']);
fprintf(fid, ['property float z', '\n']);
fprintf(fid, ['end_header', '\n']);

[h,w,~] = size(depthMap);

for i=1:h
    for j=1:w
        if(depthMap(i,j,3)~=0)
            fprintf(fid, '%f %f %f\n', depthMap(i,j,1),depthMap(i,j,2), depthMap(i,j,3));
        end
    end
end
fclose(fid);
toc
end

