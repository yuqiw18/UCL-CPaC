function distanceMatrix = DistanceMatrixRejection(rawDistanceMatrix)
    
    valueNumberToKeep = 6;
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

    % Make the matrix symmetric (two-way tracing)
    for i = 1:N
        for j = 1:N
            
            if (distanceMatrix(j,i)==0 && distanceMatrix(i,j) ~=0)
                distanceMatrix(j,i) = distanceMatrix(i,j);
                
            elseif (distanceMatrix(j,i)~=0 && distanceMatrix(i,j) ==0)
                distanceMatrix(i,j) = distanceMatrix(j,i);
            else
                
            end
        end 
    end

end