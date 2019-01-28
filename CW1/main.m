clear;
clc;

%
% cd /Users/yuqi/Desktop/UCL-CPaC/CW1/footage_raw 
% ffmpeg -i footage_%03d.png video.mp4

%
path = 'footage';
prefix = 'footage_';
first = 1; % 497
last = 657; % 657
digits = 3;
suffix = 'png';

outputPath = 'output';

verticalArtifactSequence = 2;

% footage_001.png ~ footage_657.png

% Your code here
% Image Initialisation
rawImageSequence = load_sequence(path, prefix, first, last, digits, suffix);
[vertical,horizontal,sequenceLength]=size(rawImageSequence);
imageSequence = im2double(rawImageSequence);

sceneCutFrames = [];
sceneClipPoints = [];

filter2X = [0,0,0; 0,2,0; 0,0,0];
filterLowPass = [1/9,1/9,1/9;1/9,1/9,1/9;1/9,1/9,1/9];

%% Scene Cut Detection
threshold = 45000;
tic
for i = 1+1 : sequenceLength-1
    %if (i ~= 1 && i+1<sequenceLength)       
                      
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

    %end 
end
toc

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
% Remove "redundant" cut scene frame
sceneCutFrames = sceneCutFrames(2:2:length(sceneCutFrames));

sceneClipPointsCount = length(sceneClipPoints);

% Scene-Based Operations
for p = 1 : 2 : sceneClipPointsCount-1
    for i = sceneClipPoints(p):sceneClipPoints(p+1) 
%% Global Flicker Reduction

 k = i-first+1;
% 
%     currentFrame = imageSequence(:,:,k);
%     previousFrame = imageSequence(:,:,k-1);
%     
%     intensityDiff = mean2(currentFrame) - mean2(previousFrame);  
%     previousFrame = previousFrame + intensityDiff;
%     
%     imageSequence(:,:,k-1) = previousFrame;
    
%% Blotch Correction
% S-ROD with Post-Processing

%     currentFrame = imageSequence(:,:,k);
%     previousFrame = imageSequence(:,:,k-1);
%     nextFrame = imageSequence(:,:,k+1);
%     
%     botchMask = zeros(vertical, horizontal);
%     T1 = 0.1;
%     
% for v = 2:vertical-1
%    for h = 2:horizontal-1
%        p = [previousFrame(v-1,h+1),previousFrame(v-1,h),previousFrame(v-1,h-1),nextFrame(v+1,h+1),nextFrame(v+1,h),nextFrame(v+1,h-1)];
%        if (currentFrame(v,h)-min(p)>T1)
%            disp("Test");
%        end
%    end
% end

%% Vertical Artifacts Reduction
% Medfilt with Sharpening
if (k >= sceneCutFrames(verticalArtifactSequence))
    for v = 1:vertical    
        imageSequence(v,:,k) = medfilt1(imageSequence(v,:,k),5);
    end  
end
    %imageSequence(:,:,k) = imsharpen(imageSequence(:,:,k));
    %imageSequence(:,:,k)= imfilter(imageSequence(:,:,k),filter2X) - imfilter(imageSequence(:,:,k), filterLowPass);
%% Camera Shake Calibration






        
    end 
end

% Overlay the Text
for f = 1 : length(sceneCutFrames)
    detectedFrame = imageSequence(:,:,sceneCutFrames(f));
    detectedFrame = insertText(detectedFrame,[0,0],'Scene Cut Detected','FontSize',24); 
    imageSequence(:,:,sceneCutFrames(f)) = rgb2gray(detectedFrame);
end

% Save the result
%save_sequence(imageSequence, outputPath, prefix, first, digits);
implay([im2double(rawImageSequence), imageSequence]);