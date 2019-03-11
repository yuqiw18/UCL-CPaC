% Function for computing depth map and point cloud
function [depthMap, pointCloud] = ComputeDepthMap(decodedUV, calibrationMatrix)
disp('@Computing Depth Map')
tic

    [height,width,~] = size(decodedUV);
    depthMap = zeros(height,height);
    
    % Preallocation
    pointNum = 0;
    for h = 1:height
        for w = 1:width
            if (decodedUV(h,w,1) ~= -1 && decodedUV(h,w,2) ~= -1)
                pointNum = pointNum + 1;
            end
        end
    end
    pointCloud = zeros(pointNum,3);
    
    [cameraIntrinsic,cameraExtrinsic,~,projectorIntrinsic,projectorExtrinsic,~] = LoadCalibMatrixProfile(calibrationMatrix);         
  
    % Load intrinsic matrix
    KKCamera = cameraIntrinsic;
    KKProjector = projectorIntrinsic;

    % Load rotation and translation matrix from extrinsic matrix
    RcCamera = cameraExtrinsic(1:3,1:3);
    TcCamera = cameraExtrinsic(1:3,4);
    RcProjector = projectorExtrinsic(1:3,1:3);
    TcProjector = projectorExtrinsic(1:3,4);
    
    p = 1;
    % For each point solve the linear system
    for h=1:height
        for w=1:width
            if(decodedUV(h,w,1)~=-1 && decodedUV(h,w,2)~=-1)

                % Construct A and b
                A = zeros(4,3);
                b = zeros(4,1);
                             
                pCamera = KKCamera\[w;h;1];
                pProjector = KKProjector\[decodedUV(h,w,2);decodedUV(h,w,1);1];

                A(1,:) = [RcCamera(3,1)*pCamera(1)-RcCamera(1,1), RcCamera(3,2)*pCamera(1)-RcCamera(1,2), RcCamera(3,3)*pCamera(1)-RcCamera(1,3)];
                A(2,:) = [RcCamera(3,1)*pCamera(2)-RcCamera(2,1), RcCamera(3,2)*pCamera(2)-RcCamera(2,2), RcCamera(3,3)*pCamera(2)-RcCamera(2,3)];
                A(3,:) = [RcProjector(3,1)*pProjector(1)-RcProjector(1,1), RcProjector(3,2)*pProjector(1)-RcProjector(1,2), RcProjector(3,3)*pProjector(1)-RcProjector(1,3)];
                A(4,:) = [RcProjector(3,1)*pProjector(2)-RcProjector(2,1), RcProjector(3,2)*pProjector(2)-RcProjector(2,2), RcProjector(3,3)*pProjector(2)-RcProjector(2,3)];
                b(1) = TcCamera(1)-TcCamera(3)*pCamera(1);
                b(2) = TcCamera(2)-TcCamera(3)*pCamera(2);
                b(3) = TcProjector(1)-TcProjector(3)*pProjector(1);
                b(4) = TcProjector(2)-TcProjector(3)*pProjector(2);
                
                % Solve the SLE
                x = A\b;
    
                % Acquire the Z 
                pLambda = cameraIntrinsic*cameraExtrinsic*[x;1];
                depthMap(h,w) = pLambda(3); 
                
                % Output coordinates
                pointCloud(p,:) = x;
                p = p + 1;
            else
                % Unreliable pixel
                depthMap(h,w) = -1;
            end
        end
    end

toc
end