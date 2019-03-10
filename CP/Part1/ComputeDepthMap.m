% Function for computing depth map and point cloud
function depthMap = ComputeDepthMap(decodedUV, calibrationMatrix)
disp('@Computing Depth Map')
tic

[cameraIntrinsic,cameraExtrinsic,cameraProjection,projectorIntrinsic,projectorExtrinsic,projectorProjection] = LoadCalibMatrixProfile(calibrationMatrix);
              
    [w,h,~] = size(decodedUV);
    depthMap = zeros(w,h,3);

    % Load Intrinsic Matrix
    K1 = cameraIntrinsic;
    K2 = projectorIntrinsic;

    % Compute Rotation/Translation matrices
    R1 = cameraExtrinsic(1:3,1:3);
    T1 = cameraExtrinsic(1:3,4);
    R2 = projectorExtrinsic(1:3,1:3);
    T2 = projectorExtrinsic(1:3,4);

    % For each point solve the linear system
    for i=1:w
        %disp(i/w*100);
        for j=1:h
            if(decodedUV(i,j,1)~=-1 && decodedUV(i,j,2)~=-1)
                % Compute A and b
                A = zeros(4,3);
                b = zeros(4,1);
                             
                % Convert points to normalized camera coordinates
                p1 = K1\[j;i;1];
                p2 = K2\[decodedUV(i,j,2); decodedUV(i,j,1);1];

                
                
                % Compute linear constraints -> Camera
                A(1,:) = [R1(3,1)*p1(1)-R1(1,1), R1(3,2)*p1(1)-R1(1,2), R1(3,3)*p1(1)-R1(1,3)];
                A(2,:) = [R1(3,1)*p1(2)-R1(2,1), R1(3,2)*p1(2)-R1(2,2), R1(3,3)*p1(2)-R1(2,3)];
                b(1:2) = [T1(1)-T1(3)*p1(1), T1(2)-T1(3)*p1(2)];

                % Compute linear constraints -> Projector (~ camera no.2)
                A(3,:) = [R2(3,1)*p2(1)-R2(1,1), R2(3,2)*p2(1)-R2(1,2), R2(3,3)*p2(1)-R2(1,3)];
                A(4,:) = [R2(3,1)*p2(2)-R2(2,1), R2(3,2)*p2(2)-R2(2,2), R2(3,3)*p2(2)-R2(2,3)];
                b(3:4) = [T2(1)-T2(3)*p2(1), T2(2)-T2(3)*p2(2)];

                % Compute the least squares solution
                x = A\b;

                % Store Depth
                depthMap(i,j,:) = cameraExtrinsic(1:3,1:4)*[x;1];
            end
        end
    end

toc
end