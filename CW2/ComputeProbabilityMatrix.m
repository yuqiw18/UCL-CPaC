function probabilityMatrix = ComputeProbabilityMatrix(distanceMatrix)
disp("@Compute Probability Matrix");
    tic  
    
    [h, w] = size(distanceMatrix);
    probabilityMatrix = zeros(h,w);
    
    average = 2*mean(distanceMatrix, 'all');
    
    % construct distance matrix
    for i=1:h-1
        for j=1:w-1
            probabilityMatrix(i,j) = exp(-distanceMatrix(i+1,j)/average);
        end
    end
    
    %probabilityMatrix = probabilityMatrix/sum(probabilityMatrix,'all');
    
    toc
end


