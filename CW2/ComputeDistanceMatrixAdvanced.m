% Advanced distance matrix calculation that takes both image-difference and 
% trajectory-similiarity into consideration.
function distanceMatrix = ComputeDistanceMatrixAdvanced(imageSequence, flowsData)
    
    [height, width, ~, f] = size(imageSequence);
    distanceMatrix = zeros(f,f);
    imageDifference = zeros(f,f);
    trajectorySimiliarity = zeros(f,f);
    % Construct distance matrix
    for i=1:f
        for j=1:f
            imageDifference(i,j) = sqrt(sum((imageSequence(:,:,:,i) - imageSequence(:,:,:,j)).^2,'all'));     
            if (i ~= j)
                 if (i > j)
                    k = (i-1)*(i-2)/2+j;
                    flow = flowsData(:,:,:,k);
                 else
                    k = (j-1)*(j-2)/2+i;
                    flow = -flowsData(:,:,:,k);
                 end
            else
                flow = zeros(size(flowsData(:,:,:,1)));
            end
            trajectorySimiliarity(i,j) = sqrt(sum(flow(:,:,1).^2 + flow(:,:,2).^2,'all'));          
        end
    end
    % Normalise
    imageDifference = imageDifference/max(imageDifference(:));
    trajectorySimiliarity = trajectorySimiliarity/max(trajectorySimiliarity(:));
    distanceMatrix = 2 * imageDifference + trajectorySimiliarity;
end
