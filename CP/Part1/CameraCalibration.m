warning('off','all');


%function CameraCalibration(mode)


    calibrationMode = 'synthetic';
    
    switch calibrationMode
        case 'synthetic' 
            % File Parameters
            path = 'data/synthetic_data/calibration/';
            prefix = '';
            first = 0;
            last = 5;
            digits = 4;
            suffix = 'png';
            outputPath = 'output/synthetic_reporjected_000';
            
            % Projector & Checkboard Parameters
            projectorResolution = [1024,768];
            checkboardPosition = [378,288];
            checkboardSize = 270;
            
            % Projection Matrix Estimation & Reprojection
            calibImageSequence = load_sequence_color(path, prefix, first, last, digits, suffix);   
            Reprojection(calibImageSequence, projectorResolution, checkboardPosition, checkboardSize, outputPath);
           
%             for sn=1:2:5
%             proj_name = ['data/synthetic_data/calibration/000',num2str(sn-1),'.png'];
%             cam_name = ['data/synthetic_data/calibration/000',num2str(sn),'.png'];
%             proj_img = imread(proj_name);
%             cam_img = imread(cam_name);
%             save_path = ['output/reproject_000',num2str(sn-1),'.png'];
% 
%             reprojection(proj_img, cam_img, 1024, 768, 378, 277, 270, save_path, true);
             %end
                 
        case 'real'    
        case 'capture'
        otherwise
            disp('ERROR: Undefined Calibration Mode');
    end
   
%% 1,2. Projection Matrix Estimation & Reprojection







%% 3. Point Cloud Generation & Visualisation
   
%end

function Reprojection(calibImageSequence, projectorResolution, checkboardPosition, checkboardSize, outputPath)
    
    calibImageCount = size(calibImageSequence,4);
    resWidth = projectorResolution(1);
    resHeight = projectorResolution(2);
    checkboardX = checkboardPosition(1);
    checkboardY = checkboardPosition(2);
  
    for r = 1:2:calibImageCount-1;
        figure;
        imshow(calibImageSequence(:,:,:,r));
        title('Please select 4 corners');
        [corX, corY] = getline(gcf);
        
        pts1Cart =  [corX(1:4)';corY(1:4)'];
        pts2Cart =  [checkboardX,checkboardX+checkboardSize,checkboardX+checkboardSize,checkboardX;
            checkboardY,checkboardY,checkboardY+checkboardSize,checkboardY+checkboardSize];
        
        HEst = calcBestHomography(pts1Cart, pts2Cart);
        transfrom = projective2d(HEst');
        reprojectedImage = imtransform(calibImageSequence(:,:,:,r+1),transfrom, 'XData', [1 resWidth],'YData',[1 resHeight]);
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