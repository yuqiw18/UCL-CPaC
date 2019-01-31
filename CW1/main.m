clear;
clc;

% cd /Users/yuqi/Desktop/UCL-CPaC/CW1/footage
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
pixelCount = vertical*horizontal;

imageSequence = im2double(rawImageSequence);
imageSequence2 = imageSequence;

sequenceFrameDiff = zeros(size(imageSequence));
motionMaskSequence = zeros(size(imageSequence));
blotchMaskSequence = zeros(size(imageSequence));
blotchMaskSequence2 = zeros(size(imageSequence));

sceneCutFrames = [];
sceneClipPoints = [];

sharpenFilter = [0,0,0; 0,2,0; 0,0,0];
meanFilter = [1/9,1/9,1/9;1/9,1/9,1/9;1/9,1/9,1/9];
laplacianFilter = fspecial('laplacian',0);

%% Scene Cut Detection
disp("@Scene Cut Detection");
sceneCutThreshold = 45000;
tic
for i = 2 : sequenceLength-1

        previousFrameDiff = (sum(abs(imageSequence(:,:,i)-imageSequence(:,:,i-1)),'all'));
        nextFrameDiff = (sum(abs(imageSequence(:,:,i)-imageSequence(:,:,i+1)),'all'));          
        transitionDiff = abs(previousFrameDiff - nextFrameDiff);
              
        if (transitionDiff > sceneCutThreshold)     
            sceneCutFrames = [sceneCutFrames, i+first-1];
            %disp(strcat('Scene Cut Detected @', int2str(i+first-1)));
        end
end
toc

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

% disp('Scene Cut Detected @');
% disp(sceneCutFrames);

sceneClipPointsCount = length(sceneClipPoints);

%% Global Flicker Reduction
disp("@Global Flicker Reduction");
tic
for p = 1 : 2 : sceneClipPointsCount-1
    for i = sceneClipPoints(p):sceneClipPoints(p+1)    
        k = i-first+1;    
        
        restoredFrame = imageSequence(:,:,k);
        
        [frameStart, frameEnd] = FindFrameRange(k, p, sceneClipPoints,5);       
              
%         meanHistogram = mean(imhist(imageSequence(:,:,frameStart:frameEnd)),3);
%         
%         normalizedHistogram = meanHistogram/(vertical * horizontal);
%         
%         currentFrame = histeq(currentFrame, meanHistogram);
%         
%         imageSequence(:,:,k) = currentFrame;
   
        averageIntensity = mean(imageSequence(:,:,frameStart:frameEnd),3);
        
        averageIntensity(averageIntensity>200)=200;
        
        targetHistogram = imhist(averageIntensity);
        
        restoredFrame = histeq(restoredFrame, targetHistogram);

        imageSequence(:,:,k) = restoredFrame;
        
        
%    averageIntensity = mean(imageSequence(:,:,frameStart:frameEnd),3);
    %averageIntensity = (previousFrame2 + previousFrame + currentFrame + nextFrame + nextFrame2)/5;
%    imageSequence(:,:,k) = imhistmatch(currentFrame, averageIntensity);
    end
end
toc

referenceIntensitySequence = imageSequence;

%% Blotch Correction

disp("@Blotch Correction");
tic

motionThreshold = 0.28;
blotchThreshold = 0.1;

% Motion Mask Generation
for p = 1 : 2 : sceneClipPointsCount-1
    for i = sceneClipPoints(p)+1:sceneClipPoints(p+1)
    k = i-first+1;
        sequenceFrameDiff(:,:,k) = abs(imageSequence(:,:,k) - imageSequence(:,:,k-1));   
    end
end

for p = 1 : 2 : sceneClipPointsCount-1
    for i = sceneClipPoints(p):sceneClipPoints(p+1) 
    k = i-first+1;
    [frameStart, frameEnd] = FindFrameRange(k, p, sceneClipPoints,2);
    motionMaskSequence(:,:,k) = sum(sequenceFrameDiff(:,:,frameStart:frameEnd),3);  
    end
end

% Expand the motion area so that it can cover the moving object entirely
motionMaskSequence = imfilter(motionMaskSequence, fspecial('average', 40));
motionMaskSequence(motionMaskSequence<motionThreshold )=0;
motionMaskSequence(motionMaskSequence>motionThreshold )=1;
motionMaskSequence = imdilate(motionMaskSequence, strel('disk', 24, 4));

%Blotch Correction
% blotchMaskSequence = sequenceFrameDiff;
% blotchMaskSequence(blotchMaskSequence>blotchThreshold) = 1;
% blotchMaskSequence = blotchMaskSequence .* (1-motionMaskSequence);

% Blotch Mask Generation
for p = 1 : 2 : sceneClipPointsCount-1
    for i = sceneClipPoints(p):sceneClipPoints(p+1) 
        
        k = i-first+1;

        currentFrame = imageSequence(:,:,k);
       
        if (k - 2 < sceneClipPoints(p))
            previousFrame = currentFrame;
            previousFrame2 = currentFrame;
        else
            if (k - 1 < sceneClipPoints(p))
                previousFrame = currentFrame;
                previousFrame2 = imageSequence(:,:,k-1);
            else
                previousFrame = imageSequence(:,:,k-1);
                previousFrame2 = imageSequence(:,:,k-2);
            end
        end

%         if (k + 2 > sceneClipPoints(p+1))
%             nextFrame = currentFrame;
%             nextFrame2 = currentFrame;
%         else
%             if (k + 1 > sceneClipPoints(p+1))
%                 nextFrame = currentFrame;
%                 nextFrame2 = imageSequence(:,:,k+1);
%             else
%                 nextFrame = imageSequence(:,:,k+1);
%                 nextFrame2 = imageSequence(:,:,k+2);
%             end
%         end   
        
        previousFrameDiff = abs(currentFrame - previousFrame);
        previousFrameDiff2 = abs(currentFrame - previousFrame2);

        blotchMask = zeros(size(currentFrame));
        motionMask = motionMaskSequence(:,:,k);
        
        motionMaskBin = imbinarize(motionMask);
        motionMask = imfill(motionMaskBin,'holes');    
        motionMaskSequence(:,:,k) = motionMask;

        for v = 1: vertical
            for h = 1: horizontal  
                % Mark current pixel as blotch
                if (previousFrameDiff(v,h)>blotchThreshold || previousFrameDiff2(v,h)>blotchThreshold)
                    blotchMask(v,h) = 1;
                end
                
                % Remove false blotch by checking motion mask
                if (motionMask(v,h) == 1)
                    blotchMask(v,h) = 0;
                end
            end
        end
            
        blotchMaskSequence(:,:,k)=blotchMask;

    end 
end

for p = 1 : 2 : sceneClipPointsCount-1
    for i = sceneClipPoints(p)+1:sceneClipPoints(p+1)        
        k = i-first+1;

        [frameStart, frameEnd] = FindFrameRange(k, p, sceneClipPoints,3);
        
        % Exclude current frame for blotch correction
        validFrameCount = frameEnd-frameStart;
        
        currentFrame = imageSequence(:,:,k);  
        restoredFrame = imageSequence(:,:,k);
        restoredPreviousFrame = imageSequence(:,:,k-1);

        for v = 1: vertical
            for h = 1: horizontal  
                if ( blotchMaskSequence(v,h,k) == 1)
                   % Correct blotch by getting the average value of current
                   % pixel position from previous and next few frames
                   % excluding current frame itself (blotch value = error)
                   restoredFrame(v,h) = (sum(imageSequence(v,h,frameStart:frameEnd),3) - currentFrame(v,h))/(validFrameCount);
                   
                   % Also correct previous frame because the detection
                   % compares current frames with previous frames to find
                   % blotches and previous frame may already have blotches
                   restoredPreviousFrame(v,h) = (sum(imageSequence(v,h,frameStart:frameEnd),3) - currentFrame(v,h))/(validFrameCount);
                   
                end
            end
        end
              
         imageSequence(:,:,k) = restoredFrame;
         imageSequence(:,:,k-1) = restoredPreviousFrame;
        
        % Recover the intensity gain/loss from averaging pixel values
%         imageSequence(:,:,k) = imhistmatch(restoredFrame, imageSequence(:,:,k));
%         imageSequence(:,:,k-1) = imhistmatch(restoredPreviousFrame, imageSequence(:,:,k-1));
        
    end
end


for p = 1 : 2 : sceneClipPointsCount-1
    for i = sceneClipPoints(p)+1:sceneClipPoints(p+1) 
        k = i-first+1;
        
        imageSequence(:,:,k) = imhistmatch(imageSequence(:,:,k), referenceIntensitySequence(:,:,k));
        
    end
    
end

toc

%% Vertical Artifacts Reduction
% Medfilt with Sharpening
disp("@Vertical Artifacts Reduction");
tic
for p = 1 : 2 : sceneClipPointsCount-1
    for i = sceneClipPoints(p):sceneClipPoints(p+1)    
        
        k = i-first+1;  
        
        currentFrame = imageSequence(:,:,k);       
        currentFrameFrequency = zeros(1,horizontal);
        
        if (k >= sceneCutFrames(verticalArtifactSequence))
            
            for h = 1:horizontal
                for v = 1:vertical
                    currentFrameFrequency(1,h) = currentFrameFrequency(1,h) + currentFrame(v,h);
                end
            end
            
            currentFrameFrequency = currentFrameFrequency/vertical;           
            smoothFrameFrequency = medfilt1(currentFrameFrequency,7);
            noiseFrequency = currentFrameFrequency - smoothFrameFrequency;
            
            for h = 1:horizontal
                for v = 1:vertical
                    currentFrame(v,h) = currentFrame(v,h) - noiseFrequency(1,h);
                end
            end
            
            imageSequence(:,:,k) = currentFrame;

%             for v = 1:vertical
%                 for m = 5:-1:1
%                 imageSequence(v,:,k) = medfilt1(imageSequence(v,:,k),m);
%                 end
%             end

            % Recover features from blurring
            % Laplacian
%             imageSequenceEdge(:,:,k) = imfilter(imageSequence(:,:,k), laplacianFilter, 'replicate');
%             imageSequence(:,:,k) = imageSequence(:,:,k) - imageSequenceEdge(:,:,k);

            % Sharpening
%             imageSequence(:,:,k) = imsharpen(imageSequence(:,:,k));
%             imageSequence(:,:,k)= imfilter(imageSequence(:,:,k),sharpenFilter) - imfilter(imageSequence(:,:,k), meanFilter);
        end
    end
end
toc

%% Camera Shake Calibration

%         currentFrame = uint8(currentFrame);
%         previousFrame = uint8(previousFrame);
%         
%         % detect corners of prev and curr frame
%         currentPoints = detectFASTFeatures(currentFrame, 'MinContrast', 0.1);
%         previousPoints = detectFASTFeatures(previousFrame, 'MinContrast', 0.1);
%         
%         % extract descriptors from corners
%         [currentFeatures, currentPoints] = extractFeatures(currentFrame, currentPoints);
%         [previousFeatures, previousPoints] = extractFeatures(previousFrame, previousPoints);
%         
%         % match features
%         indexPairs = matchFeatures(currentFeatures, previousFeatures);
%         currentPoints = currentPoints(indexPairs(:, 1), :);
%         previousPoints = previousPoints(indexPairs(:, 2), :);
%         
%         if size(currentPoints, 1) >= 3
%             tform = estimateGeometricTransform(previousPoints, currentPoints, 'affine');
%             imageSequence(:,:,k) = imwarp(previousFrame, tform, 'OutputView',imref2d(size(currentFrame)));
%         end

%%
% Overlay the Text
for f = 1 : length(sceneCutFrames)
    detectedFrame = imageSequence(:,:,sceneCutFrames(f));
    detectedFrame = insertText(detectedFrame,[0,0],'Scene Cut Detected','FontSize',24); 
    imageSequence(:,:,sceneCutFrames(f)) = rgb2gray(detectedFrame);
end

% Save the result
%save_sequence(imageSequence, outputPath, prefix, first, digits);

% Frame by frame comparasion
implay([blotchMaskSequence,motionMaskSequence,sequenceFrameDiff, im2double(rawImageSequence), imageSequence]);

%%
% Functions
% Find closest start and end frames
function [frameStart, frameEnd] = FindFrameRange(currentFrame, currentSequenceNumber, sceneClipPoints, n)
    
startLimit = sceneClipPoints(currentSequenceNumber);
    endLimit = sceneClipPoints(currentSequenceNumber+1);
    
    for i = n:-1:0
        if (currentFrame - i < startLimit)
            
        else
            frameStart = currentFrame - i;
            break;
        end
    end
    
    for i = n:-1:0
        if (currentFrame + i > endLimit)
            
        else
            frameEnd = currentFrame + i;
            break;
        end
    end
end