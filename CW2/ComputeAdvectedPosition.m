% Function for calculating advected position using optical flow
% The algorithm is based on coursework instruction
function advectedPositionValue = ComputeAdvectedPosition(startPoint, paths, flowsData)

    pathCount = length(paths);
    advectedPositionValue = zeros(pathCount,2);
    
    % For each path
    for i = 1:pathCount
        
        currentPath = paths{i};  
        currentLocation = startPoint;
        
        currentPosition = round(currentLocation);      
        currentPositionX = currentPosition(1);         
        currentPositionY = currentPosition(2);
        
        % For each node on this path
        if length(currentPath)>1
            for f=1:length(currentPath)-1
                
                currentFrameIndex = currentPath(f);       
                nextFrameIndex = currentPath(f+1);
                
                if (currentFrameIndex>nextFrameIndex)
                    k = (currentFrameIndex-1)*(currentFrameIndex-2)/2 + nextFrameIndex;    
                    flowX = flowsData(currentPositionY,currentPositionX,1,k);
                    flowY = flowsData(currentPositionY,currentPositionX,2,k);
                elseif (currentFrameIndex<nextFrameIndex)
                    k = (nextFrameIndex-1)*(nextFrameIndex-2)/2 + currentFrameIndex;
                    flowX = -flowsData(currentPositionY,currentPositionX,1,k);
                    flowY = -flowsData(currentPositionY,currentPositionX,2,k);
                else
                    flowX = 0;
                    flowY = 0;
                end
                
                % Accumulate the position value
                currentLocation = currentLocation + [flowX, flowY];  
            end   
            advectedPositionValue(i,:) = currentLocation;
        else  
            advectedPositionValue(i,:) = currentLocation;          
        end    
    end 
end
