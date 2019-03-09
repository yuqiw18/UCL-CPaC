%function CameraCalibration(mode)

    calibrationMode = 'synthetic';
    
    switch calibrationMode
        case 'synthetic' 
%             % File Parameters
%             path = 'data/synthetic_data/calibration/';
%             prefix = '';
%             first = 0;
%             last = 5;
%             digits = 4;
%             suffix = 'png';
%             outputPath = 'output/synthetic_reporjected_000';
%             
%             % Projector & Checkboard Parameters
%             projectorResolution = [1024,768];
%             checkboardPosition = [378,288];
%             checkboardSize = 270;
%             
%             % Projection Matrix Estimation & Reprojection
%             calibImageSequence = double(load_sequence(path, prefix, first, last, digits, suffix));   
%             Reprojection2(calibImageSequence, projectorResolution, checkboardPosition, checkboardSize, outputPath);
            for sn=1:2:5
            proj_name = ['data/synthetic_data/calibration/000',num2str(sn-1),'.png'];
            cam_name = ['data/synthetic_data/calibration/000',num2str(sn),'.png'];
            proj_img = imread(proj_name);
            cam_img = imread(cam_name);
            save_path = ['output/reproject_000',num2str(sn-1),'.png'];

            reprojection(proj_img, cam_img, 1024, 768, 378, 277, 270, save_path, true);
%             imwrite(cam_img_out, save_path);
             end
                 
        case 'real'
            for rn=1:9
            proj_name = ['data/real_calibration/IMG_932',num2str(rn),'.jpg'];
            cam_name = ['data/real_calibration/IMG_932',num2str(rn),'.jpg'];
            proj_img = imread(proj_name);
            cam_img = imread(cam_name);
            save_path = ['data/real_calibration/reproject_932',num2str(rn),'.jpg'];
            
            [proj_img_out, cam_img_out] = reprojection(proj_img, cam_img, 1024, 768, 518, 120, 299, true);
            imwrite(cam_img_out,save_path);
            end
        case 'capture'
            for on=1:6
            proj_name = ['data/own_calibration/IMG_28',num2str(on+87),'.jpg'];
            cam_name = ['data/own_calibration/IMG_28',num2str(on+87),'.jpg'];
            proj_img = imread(proj_name);
            cam_img = imread(cam_name);
            save_path = ['data/own_calibration/reproject_28',num2str(on+87),'.jpg'];
            
            [proj_img_out, cam_img_out] = reprojection(proj_img, cam_img, 1920, 1080, 705, 285, 512, true);
            imwrite(cam_img_out,save_path);
            end
        otherwise
            disp('ERROR: Undefined Calibration Mode');
    end
   
%% 1,2. Projection Matrix Estimation & Reprojection







%% 3. Point Cloud Generation & Visualisation
   
%end

function Reprojection2(calibImageSequence, projectorResolution, checkboardPosition, checkboardSize, outputPath)
    
    calibImageCount = size(calibImageSequence,3);
    resWidth = projectorResolution(1);
    resHeight = projectorResolution(2);
    checkboardX = checkboardPosition(1);
    checkboardY = checkboardPosition(2);
  
    for r = 1:2:calibImageCount-1;
        figure;
        imshow(uint8(calibImageSequence(:,:,r)));
        title('Please select 4 corners');
        [corX, corY] = getline(gcf);
        
        pts1Cart =  [corX(1:4)';corY(1:4)'];
        pts2Cart =  [checkboardX,checkboardX+checkboardSize,checkboardX+checkboardSize,checkboardX;
            checkboardY,checkboardY,checkboardY+checkboardSize,checkboardY+checkboardSize];
        
        HEst = calcBestHomography(pts1Cart, pts2Cart);
        transfrom = projective2d(HEst');
        reprojectedImage = imwarp(calibImageSequence(:,:,r+1),transfrom);
        imwrite(reprojectedImage,strcat(outputPath, num2str(r+1), '.png'));
        close all;
    end
    
end



function reprojection(I_0, I_1, w_p, h_p, x_0, y_0, dx, path, resize)

% =========== Parameters ============
% I_0 : Input Image (projected)
% I_1 : Output Image (to be rectified)
% [w_p,h_p] : Output Image Size
% [x_0,y_0] : Checkerboard Offset
% dx : Checkerboard Size (pixels)
% ====================================

% Enter Four Corners
figure
imshow(I_0);
title('Select Corners');
[x, y] = getline(gcf);
close all;

% Compute Homography
pts1Cart =  [x(1:4)';y(1:4)'];
pts2Cart =  [x_0,  x_0+dx, x_0+dx, x_0;...
             y_0,  y_0,    y_0+dx, y_0+dx];

% Estimate Homography
HEst = calcBestHomography(pts1Cart, pts2Cart);

% Reproject Image
t = maketform('projective',HEst');
if(resize)
    I_0_out = imtransform(I_0,t,'XData',[1 w_p], 'YData', [1 h_p]);
    I_1_out = imtransform(I_1,t,'XData',[1 w_p], 'YData', [1 h_p]);
else
    I_0_out = imtransform(I_0,t);
    I_1_out = imtransform(I_1,t);
end

% Save Image
imwrite(I_1_out,path);
end

%% Plotting
% figure
% subplot(2,2,1)
% imshow(I_0); title('Original Projected Image');
% subplot(2,2,2)
% imshow(I_1); title('Original Checkerboard Image');
% subplot(2,2,3)
% imshow(I_0_out); title('Rectified Projected Image');
% subplot(2,2,4)
% imshow(I_1_out); title('Reprojected Checkerboard Image');


%% ==========================================================================
function H = calcBestHomography(pts1Cart, pts2Cart)

pts1Cart = [pts1Cart; ones(1,size(pts1Cart,2))];
pts2Cart = [pts2Cart; ones(1,size(pts2Cart,2))];

A = zeros(2*size(pts1Cart,2),9);
for i = 1:size(pts1Cart,2)
    ui = pts1Cart(1,i);
    vi = pts1Cart(2,i);
    xi = pts2Cart(1,i);
    yi = pts2Cart(2,i);
    A(2*i-1,:) = [0,0,0,-ui,-vi,-1,yi*ui,yi*vi,yi];
    A(2*i,:) = [ui,vi,1,0,0,0,-xi*ui,-xi*vi,-xi];
end

h = solveAXEqualsZero(A);
H = reshape(h,[3,3])';
end
%% ==========================================================================
function x = solveAXEqualsZero(A)
[~,~,V] = svd(A);
x = V(:,end);
end