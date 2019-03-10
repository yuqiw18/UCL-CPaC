warning('off','all');
clear;

%function CameraCalibration(mode)
calibrationMode = 'real';

switch calibrationMode
    case 'synthetic' 
        % File Parameters
        path = 'data/synthetic_data/calibration/';
        prefix = '';
        first = 0;
        last = 5;
        digits = 4;
        suffix = 'png';
        outputPath = 'output/synthetic_reprojected_';

        % Projector & Checkboard Parameters
        projectorResolution = [1024,768];
        checkboardPosition = [378,277];
        checkboardSize = 270;

        % Projection Matrix Estimation & Reprojection
        imageSequence = load_sequence_color(path, prefix, first, last, digits, suffix);   
        calibImageSequence = imageSequence(:,:,:,1:2:5);
        targetImageSequence = imageSequence(:,:,:,2:2:6);
        clear imageSequence;
        
        Reprojection(calibImageSequence, targetImageSequence, projectorResolution, checkboardPosition, checkboardSize, outputPath);

    case 'real'    
        % File Parameters
        path = 'data/real_data/real_calibration/';
        prefix = 'IMG_';
        first = 9321;
        last = 9329;
        digits = 4;
        suffix = 'jpg';
        outputPath = 'output/real_reprojected_';

        % Projector & Checkboard Parameters
        projectorResolution = [1024,768];
        checkboardPosition = [518,120];
        checkboardSize = 299;

        % Projection Matrix Estimation & Reprojection
        calibImageSequence = load_sequence_color(path, prefix, first, last, digits, suffix);
        targetImageSequence = calibImageSequence;
        Reprojection(calibImageSequence, targetImageSequence, projectorResolution, checkboardPosition, checkboardSize, outputPath);
        
    case 'capture'
    otherwise
        disp('ERROR: Undefined Calibration Mode');
end
  
function Reprojection(calibImageSequence,targetImageSequence, projectorResolution, checkerboardPosition, checkerboardSize, outputPath)
    
    calibImageCount = size(calibImageSequence,4);
    resWidth = projectorResolution(1);
    resHeight = projectorResolution(2);
    checkerboardX = checkerboardPosition(1);
    checkerboardY = checkerboardPosition(2);
    checkerboardCorner = [checkerboardX, checkerboardX+checkerboardSize, checkerboardX+checkerboardSize, checkerboardX;
                               checkerboardY, checkerboardY, checkerboardY+checkerboardSize, checkerboardY+checkerboardSize];
                           
    for r = 1:calibImageCount
        figure;
        imshow(calibImageSequence(:,:,:,r));
        title('Please select 4 corners');
        [corX, corY] = getline(gcf);
        selectedCorner = [corX, corY];
        
        H = EstimateHomography(checkerboardCorner, selectedCorner');
        transfrom = maketform('projective',H');
        reprojectedImage = imtransform(im2double(targetImageSequence(:,:,:,r)),transfrom, 'XData', [1, resWidth],'YData',[1, resHeight]);
        imwrite(reprojectedImage,strcat(outputPath, num2str(r), '.jpg'));
        close all;
    end
    
end