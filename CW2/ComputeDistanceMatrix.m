% Function for computing the distance matrix (basic)
% Basic distance matrix only takes the image difference into account.
function distanceMatrix = ComputeDistanceMatrix(imageSequence)
disp("@Computing Distance Matrix");
tic  
    [~, ~, ~, f] = size(imageSequence);
    distanceMatrix = zeros(f, f);
    
    % Construct Distance Matrix
    for i=1:f
        for j=1:f
            currentDistance = sqrt(sum((imageSequence(:,:,:,i) - imageSequence(:,:,:,j)).^2,'all'));
            distanceMatrix(i,j) = currentDistance;
        end
    end
    
    % Normalise
    distanceMatrix = distanceMatrix/(max(distanceMatrix(:)));
toc
end





