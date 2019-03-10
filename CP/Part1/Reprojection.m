
function Reprojection(calibImageSequence,targetImageSequence, projectorResolution, checkerboardPosition, checkerboardSize, outputPath)
    
    calibImageCount = size(calibImageSequence,4);
    resWidth = projectorResolution(1);
    resHeight = projectorResolution(2);
    checkerboardX = checkerboardPosition(1);
    checkerboardY = checkerboardPosition(2);
    
    % Define checkerboard corners (clockwise-> top-left, top-right, bottom-right, bottom-left)
    checkerboardCorner = [checkerboardX, checkerboardX+checkerboardSize, checkerboardX+checkerboardSize, checkerboardX;
                               checkerboardY, checkerboardY, checkerboardY+checkerboardSize, checkerboardY+checkerboardSize];
                           
    for r = 1:calibImageCount
        
        % Select 4 corners (clockwise-> top-left, top-right, bottom-right, bottom-left)
        figure;
        imshow(calibImageSequence(:,:,:,r));
        title('Please select 4 corners');
        [corX, corY] = getline(gcf);
        selectedCorner = [corX, corY];
        
        % Estimate the homography from selected corners
        H = EstimateHomography(checkerboardCorner, selectedCorner');
        
        % Reproject with homography
        transform = projective2d(H');
        reprojectedImage = imwarp(im2double(targetImageSequence(:,:,:,r)), transform, 'OutputView', imref2d([resHeight,resWidth]));
        
        % Save the result for next calibration steps
        imwrite(reprojectedImage,strcat(outputPath, num2str(r), '.jpg'));
        close all;
    end
    
end