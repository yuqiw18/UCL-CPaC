function [closestIndex, closestX, closestY] = FindShortestAdvectedPath(startPoint, endPoint, paths, flowsData)

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
                
                % The width and height of optical flow is actually the
                % height and width of the image data. WidthOF = HeightI,
                % HeightOF = WidthI
                % Thus we need to swap X and Y coordinate here
                if (currentFrameIndex>nextFrameIndex)
                    k = (currentFrameIndex-1)*(currentFrameIndex-2)/2 + nextFrameIndex;    
                    flowX = flowsData(currentPositionY,currentPositionX,1,k);
                    flowY = flowsData(currentPositionY,currentPositionX,2,k);
                    
                elseif (currentFrameIndex<nextFrameIndex)
                    
                    k = (nextFrameIndex-1)*(nextFrameIndex-2)/2+currentFrameIndex;
                    flowX = -flowsData(currentPositionY,currentPositionX,1,k);
                    flowY = -flowsData(currentPositionY,currentPositionX,2,k);
                    
                else
                    flowX = 0;
                    flowY = 0;
                end
                
                currentLocation = currentLocation + [flowX, flowY];
            end   
            
            %
            positionFlowValue(i,:) = currentLocation;
            
        else 
            
            %
            positionFlowValue(i,:) = currentLocation;          
        end    
    end
    
    % Find closest advected position
    advectedDistances = zeros(pathCount,1);
    for i = 1:pathCount 
        advectedDistances(i) = sqrt(sum((positionFlowValue(i,:)-endPoint).^2,'all'));
    end
    
    [~, closestIndex] = min(advectedDistances);
    
    closestX = positionFlowValue(closestIndex,1);
    closestY = positionFlowValue(closestIndex,2);
    
    disp(closestIndex);
end
