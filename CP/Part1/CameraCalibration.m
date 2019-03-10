warning('off','all');
clear all;

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
        outputPath = 'output/synthetic_reprojected_000';

        % Projector & Checkboard Parameters
        projectorResolution = [1024,768];
        checkboardPosition = [378,277];
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
function Reprojection(calibImageSequence, projectorResolution, checkerboardPosition, checkerboardSize, outputPath)
    
    calibImageCount = size(calibImageSequence,4);
    resWidth = projectorResolution(1);
    resHeight = projectorResolution(2);
    checkerboardX = checkerboardPosition(1);
    checkerboardY = checkerboardPosition(2);
    checkerboardCorner = [checkerboardX, checkerboardX+checkerboardSize, checkerboardX+checkerboardSize, checkerboardX;
                               checkerboardY, checkerboardY, checkerboardY+checkerboardSize, checkerboardY+checkerboardSize];
                           
    for r = 1:2:calibImageCount-1
        figure;
        imshow(calibImageSequence(:,:,:,r));
        title('Please select 4 corners');
        [corX, corY] = getline(gcf);
        selectedCorner = [corX, corY];
        
        H = EstimateHomography(checkerboardCorner, selectedCorner');
        %H = get_homography(selectedCorner', checkerboardCorner);
        transfrom = maketform('projective',H');
        reprojectedImage = imtransform(im2double(calibImageSequence(:,:,:,r+1)),transfrom, 'XData', [1, resWidth],'YData',[1, resHeight]);
        imwrite(reprojectedImage,strcat(outputPath, num2str(r+1), '.jpg'));
        close all;
    end
    
end