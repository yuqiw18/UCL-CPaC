% Function for generating image from given paths (no duplication between 2
% frames) and index list as well
function [outputIndex, outputImageSequence] = ConvertPathsToImageSequence(paths, sourceImageSequence)
disp("@Converting Paths to Image Sequence");
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
    [height,width,channel,~] = size(sourceImageSequence);
    outputImageSequence = zeros(height, width, channel, length(outputIndex));
    for f = 1:length(outputIndex)
       outputImageSequence(:,:,:,f) = sourceImageSequence(:,:,:,outputIndex(f)); 
    end  
toc
end