% Generate image sequence from given paths (no duplication between 2 frames)
function outputImageSequence = ConvertPathsToImageSequence(paths, sourceImageSequence)
disp("@Convert Paths to Image Sequence");
tic  
    % Convert paths to index
    outputIndex = [];
    lastIndex = 0;
    for p = 1:length(paths)
        path = paths{p};
        for i = 1:length(path)
            if (path(i) ~= lastIndex)
                outputIndex = [outputIndex path(i)];
                lastIndex = path(i);
            end
        end
    end
  
    % Extract frames using index
    [h,w,c,~] = size(sourceImageSequence);
    
    outputImageSequence = zeros(h, w, c, length(outputIndex));
    for f = 1:length(outputIndex)
       outputImageSequence(:,:,:,f) = sourceImageSequence(:,:,:,outputIndex(f)); 
    end  
toc
end