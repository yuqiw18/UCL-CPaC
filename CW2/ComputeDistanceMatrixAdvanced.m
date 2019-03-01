% Advanced distance matrix calculation that takes both image-difference and 
% trajectory-similiarity into consideration.
function distanceMatrix = ComputeDistanceMatrixAdvanced(imageSequence, flowsData)
    
    [height, width, ~, f] = size(imageSequence);
    distanceMatrix = zeros(f,f);
    % Construct distance matrix
    for i=1:f
        for j=1:f
            currentDistance = (imageSequence(:,:,:,j) - imageSequence(:,:,:,i)).^2;
            currentDistance = sum(currentDistance(:))/(height*width);
            
            if (i ~= j)
                 if i>j
                    k = (i-1)*(i-2)/2+j;
                    flowX = flowsData(:,:,1,k).^2;
                    flowY = flowsData(:,:,2,k).^2;
                 else
                    k = (j-1)*(j-2)/2+i;
                    flowX = flowsData(:,:,1,k).^2;
                    flowY = flowsData(:,:,2,k).^2;
                 end
                 flowXDiff = sum(flowX(:))/(height*width);
                 flowYDiff = sum(flowY(:))/(height*width);
            else
                flowXDiff = 0;
                flowYDiff = 0;
            end
            distanceMatrix(i,j) = sqrt(currentDistance+flowXDiff+flowYDiff);
        end
    end
    
    % Normalise
    distanceMatrix = distanceMatrix/max(distanceMatrix(:));
end
