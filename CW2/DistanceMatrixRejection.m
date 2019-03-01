function distanceMatrix = DistanceMatrixRejection(rawDistanceMatrix)
    
    valueNumberToKeep = 8;
    [N,~] = size(rawDistanceMatrix);
    
    distanceMatrix = zeros(size(rawDistanceMatrix));

    for i = 1:N
        
        currentRow = rawDistanceMatrix(i,1:N);
        valueList = currentRow(currentRow>0);
  
        sortedValue = sort(valueList);
        
        for v = 1:valueNumberToKeep
            index = find(currentRow == sortedValue(v));
            distanceMatrix(i,index) = sortedValue(v);
        end
        
    end

    for i = 1:N
        for j = 1:N
            distanceMatrix(j,i) = distanceMatrix(i,j);
        end 
    end
    
    imshow(distanceMatrix);
end