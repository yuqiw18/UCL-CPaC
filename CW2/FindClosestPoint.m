% Function for finding the closest point
% Pick the path in Paths whose advected location comes closest to the
% end point
function [closestIndex, closestPoint] = FindClosestPoint(endPoint, advectedPositionValues)

    valueCount = length(advectedPositionValues);
    
    % Compute the distance between end point and each advected location 
    advectedDistances = zeros(valueCount,1);
    for i = 1:valueCount 
        advectedDistances(i) = sqrt(sum((advectedPositionValues(i,:)-endPoint).^2,'all'));
    end
    
    % Find the index and coordinates of the minimum distance
    [~, closestIndex] = min(advectedDistances);
    closestPoint = [advectedPositionValues(closestIndex,1), advectedPositionValues(closestIndex,2)];
    
 end