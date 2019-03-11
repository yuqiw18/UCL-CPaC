warning('off','all');
clear;

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
        % File Parameters
        path = 'data/capture_data/calibration/';
        prefix = 'IMG_';
        first = 3247;
        last = 3255;
        digits = 4;
        suffix = 'jpg';
        outputPath = 'output/capture_reprojected_';

        % Projector & Checkboard Parameters
        projectorResolution = [864,539];
        checkboardPosition = [434,236];
        checkboardSize = 232;

        % Projection Matrix Estimation & Reprojection
        calibImageSequence = load_sequence_color(path, prefix, first, last, digits, suffix);
        targetImageSequence = calibImageSequence;
        Reprojection(calibImageSequence, targetImageSequence, projectorResolution, checkboardPosition, checkboardSize, outputPath);
    otherwise
        disp('ERROR: Undefined Calibration Mode');
end
