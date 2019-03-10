function H = EstimateHomography(checkerboardCorner, selectedCorner)

    [~,num] = size(selectedCorner);
    % Convert the cordinates into Cartesian system
    pts1Cart = [selectedCorner; ones(1,num)];
    pts2Cart = [checkerboardCorner; ones(1,num)];

    A = zeros(2*num,9);

    for i = 1:num
        ui = pts1Cart(1,i);
        vi = pts1Cart(2,i);
        xi = pts2Cart(1,i);
        yi = pts2Cart(2,i);
        A(2*i-1,:) = [0,0,0,-ui,-vi,-1,yi*ui,yi*vi,yi];
        A(2*i,:) = [ui,vi,1,0,0,0,-xi*ui,-xi*vi,-xi];
    end

    [~,~,V] = svd(A);
    h = V(:,size(V,2));
    H = reshape(h,[3,3])';
    
end