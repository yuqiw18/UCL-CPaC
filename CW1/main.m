clear;
clc;

%
% cd /Users/yuqi/Desktop/UCL-CPaC/CW1/footage_raw 
% ffmpeg -i footage_%03d.png video.mp4

%
path = 'footage';
prefix = 'footage_';
first = 497;
last = 657;
digits = 3;
suffix = 'png';

outputPath = 'output';

% footage_001.png ~ footage_657.png

% Your code here
%% Image Initialisation
rawImageSequence = load_sequence(path, prefix, first, last, digits, suffix);
[~,~,sequenceLength]=size(rawImageSequence);
imageSequence = im2double(rawImageSequence);

sceneCutFrames = [];
sceneClipPoints = [];

%% Scene Cut Detection
threshold = 45000;
for i = 1 : sequenceLength 
    if (i ~= 1 && i+1<sequenceLength)       
                      
        %
        previousFrameDiff = (sum(abs(imageSequence(:,:,i)-imageSequence(:,:,i-1)),'all'));
        nextFrameDiff = (sum(abs(imageSequence(:,:,i)-imageSequence(:,:,i+1)),'all'));          
        transitionDiff = abs(previousFrameDiff - nextFrameDiff);
              
        if (transitionDiff > threshold)     
            sceneCutFrames = [sceneCutFrames, i+first-1];
            %disp(strcat('Scene Cut Detected @', int2str(i+first-1)));
        end
        
        % Debug
%         disp(strcat("Previous Frame Diff: ",num2str(abs(previousFrameDiff))));
%         disp(strcat("Next Frame Diff: ",num2str(abs(nextFrameDiff))));
%         disp(strcat("Transition Diff: ", num2str(transitionDiff)));
%         disp(" ");        
    end 
end

 disp('Scene Cut Detected @');
 disp(sceneCutFrames);

% Footages between 2 points will be a continuous scene e.g. 
% 1 - 1xx is one scene
% 1xx + 1 - 2xx is another scene

sceneClipPoints = sceneCutFrames;
% Append first frame
sceneClipPoints = [first, sceneClipPoints];
% Append last frame
sceneClipPoints = [sceneClipPoints, last];

%% Scene-Based Operations
for p = 1 : 2 : length(sceneClipPoints)-1
    for i = sceneClipPoints(p)+1:sceneClipPoints(p+1)    
% Global Flicker Correction

k = i-first+1;


    currentFrame = imageSequence(:,:,k);
    previousFrame = imageSequence(:,:,k-1);
    
    intensityDiff = mean2(currentFrame) - mean2(previousFrame);  
    previousFrame = previousFrame + intensityDiff;
    
    imageSequence(:,:,k-1) = previousFrame;
    
% Blotch Correction






% Vertical Artifacts Correction



imageSequence(:,:,k) = medfilt1(imageSequence(:,:,k), 5);


% Camera Shake Calibration



        
        
        
    end 
end
%% Save the result
save_sequence(imageSequence, outputPath, prefix, first, digits);