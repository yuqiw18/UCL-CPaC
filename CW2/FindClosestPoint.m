function [closestIndex, closestPoint] = FindClosestPoint(endPoint, advectedPositionValues)

    valueCount = length(advectedPositionValues);
    
    advectedDistances = zeros(valueCount,1);
    for i = 1:valueCount 
        advectedDistances(i) = sqrt(sum((advectedPositionValues(i,:)-endPoint).^2,'all'));
    end
    
    [~, closestIndex] = min(advectedDistances);
    closestPoint = [advectedPositionValues(closestIndex,1), advectedPositionValues(closestIndex,2)];
    
 end