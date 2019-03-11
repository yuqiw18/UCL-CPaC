% Function for estimate the homography
function H = EstimateHomography(checkerboardCorner, selectedCorner)

    % Convert to homogeneous image coordinates
    % (x,y) -> [x y 1];
    hSelectedCorner = [selectedCorner; ones(1,4)];
    hCheckerboardCorner = [checkerboardCorner; ones(1,4)];

    % Construct A
    % 4 paris of correspondance points
    A = zeros(8,9);

    for n = 1:4 
        xS = hSelectedCorner(1,n);
        yS = hSelectedCorner(2,n);
        
        xC = hCheckerboardCorner(1,n);
        yC = hCheckerboardCorner(2,n);
        
        A(2*n-1,:) = [0,0,0,-xS,-yS,-1,yC*xS,yC*yS,yC];
        A(2*n,:) = [-xS,-yS,-1,0,0,0,xC*xS,xC*yS,xC];
    end
    
    % Solve the SVD where AH = 0 and A=USV';
    % Select the last singular vector of V as the solution to H
    [~,~,V] = svd(A); 
    s = V(:,size(V,2))';
    
    % Reshape from a 1x9 matrix to 3x3 matrix for use in transformation
    H = reshape(s,[3,3]);
    
end