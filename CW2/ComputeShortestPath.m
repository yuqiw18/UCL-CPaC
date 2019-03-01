function [closestIndex,closestX, closestY] = ComputeShortestPath(startPoint, endPoint, paths, flowsData)

    pathCount = length(paths);
    positionFlowValue = zeros(pathCount,2);
    
    for i = 1:pathCount
        
        currentPath = paths{i};   
        currentLocation = startPoint;
        
        if length(currentPath)>1
            
            for f=1:length(currentPath)-1
                
                currentFrameIndex = currentPath(f);
                nextFrameIndex = currentPath(f+1);
                
                currentPosition = round(currentLocation); 
                
                currentPositionX = currentPosition(1);
                currentPositionY = currentPosition(2);
                 
                if (currentFrameIndex>nextFrameIndex)
                    k = (currentFrameIndex-1)*(currentFrameIndex-2)/2 + nextFrameIndex;    
                    flowX = flowsData(currentPositionX,currentPositionY,1,k);
                    flowY = flowsData(currentPositionX,currentPositionY,2,k);
                    
                elseif (currentFrameIndex<nextFrameIndex)
                    
                    k = (nextFrameIndex-1)*(nextFrameIndex-2)/2+currentFrameIndex;
                    flowX = -flowsData(currentPositionX,currentPositionY,1,k);
                    flowY = -flowsData(currentPositionX,currentPositionY,2,k);
                    
                else
                    flowX = 0;
                    flowY = 0;
                end
                currentLocation = currentLocation + [flowY, flowX];
            end   
            
            %
            positionFlowValue(i,:) = currentLocation;
            
        else 
            
            %
            positionFlowValue(i,:) = currentLocation;          
        end    
    end
    
    advectedDistance = sqrt(sum((positionFlowValue-repmat(endPoint,length(paths),1)).^2,2)); 
    
    [~, closestIndex] = min(advectedDistance);
    
    closestX = positionFlowValue(closestIndex,2);
    closestY = positionFlowValue(closestIndex,1);
end
