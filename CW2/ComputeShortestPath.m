function [closestIndex,closestX, closestY] = ComputeShortestPath(startPoint, endPoint, paths, flowsData)

    last_pixs = zeros(length(paths),2);
    
    for i = 1:length(paths)
        
        currentPath = paths{i};
        currentLocation = startPoint;
        
        if length(currentPath)>1
            
            for f=1:length(currentPath)-1
                
                currentFrameIndex = currentPath(f);
                nextFrameIndex = currentPath(f+1);
                
                current_pos = round(currentLocation);    
                 
                if (currentFrameIndex>nextFrameIndex)
                    k = (currentFrameIndex-1)*(currentFrameIndex-2)/2 + nextFrameIndex;    
                    flowX = flowsData(current_pos(1),current_pos(2),1,k);
                    flowY = flowsData(current_pos(1),current_pos(2),2,k);
                else
                    k = (nextFrameIndex-1)*(nextFrameIndex-2)/2+currentFrameIndex;
                    flowX = -flowsData(current_pos(1),current_pos(2),1,k);
                    flowY = -flowsData(current_pos(1),current_pos(2),2,k);
                end
                currentLocation = currentLocation + [flowY, flowX];
            end   
            last_pixs(i,:) = currentLocation;
        else
            last_pixs(i,:) = currentLocation;
        end    
    end
    
    dist_toExpect = sqrt(sum((last_pixs-repmat(endPoint,length(paths),1)).^2,2)); 
    
    [min_dist, closestIndex] = min(dist_toExpect);
    
    closestX = last_pixs(closestIndex,2);
    closestY = last_pixs(closestIndex,1);
end
