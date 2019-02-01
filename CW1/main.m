clear;
clc;

path = 'footage';
prefix = 'footage_';
first = 1; % 497
last = 498; % 657
digits = 3;
suffix = 'png';

outputPath = 'output';

verticalArtifactSequence = 2;

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

sceneCutFrames = [];
sceneClipPoints = [];

sharpenFilter = [0,0,0; 0,2,0; 0,0,0];
meanFilter = [1/9,1/9,1/9;1/9,1/9,1/9;1/9,1/9,1/9];
laplacianFilter = fspecial('laplacian',0);

%% Scene Cut Detection
disp("@Scene Cut Detection");

sceneCutThreshold = 45000/pixelCount;
tic
for i = 2 : sequenceLength-1

        previousFrameDiff = (sum(abs(imageSequence(:,:,i)-imageSequence(:,:,i-1)),'all'));
        nextFrameDiff = (sum(abs(imageSequence(:,:,i)-imageSequence(:,:,i+1)),'all'));          
        transitionDiff = abs(previousFrameDiff - nextFrameDiff)/pixelCount;
              
        if (transitionDiff > sceneCutThreshold)     
            sceneCutFrames = [sceneCutFrames, i+first-1];
            %disp(strcat(' - Scene Cut Detected @', int2str(i+first-1)));
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

        averageIntensity = mean(imageSequence(:,:,frameStart:frameEnd),3);
        
        averageIntensity(averageIntensity>200)=200;
        
        targetHistogram = imhist(averageIntensity);
        
        restoredFrame = histeq(restoredFrame, targetHistogram);

        imageSequence(:,:,k) = restoredFrame;
        
    end
end
toc

%% Blotch Correction

disp("@Blotch Correction");
tic

motionThreshold = 0.28;
motionThresholdRefine = 0.9;
blotchThreshold = 0.1;

% Motion Mask Generation
for p = 1 : 2 : sceneClipPointsCount-1
    for i = sceneClipPoints(p)+1:sceneClipPoints(p+1)
    k = i-first+1;
        sequenceFrameDiff(:,:,k) = abs(imageSequence(:,:,k) - imageSequence(:,:,k-1));   
    end
    
    for i = sceneClipPoints(p):sceneClipPoints(p+1) 
    k = i-first+1;
    [frameStart, frameEnd] = FindFrameRange(k, p, sceneClipPoints,2);
    motionMaskSequence(:,:,k) = sum(sequenceFrameDiff(:,:,frameStart:frameEnd),3);  
    end  
end

expansionFilter = fspecial('average', 35);
dilateStructure = strel('disk', 28, 4);

% Expand the motion area so that it can cover the moving object entirely
motionMaskSequenceRaw = imfilter(motionMaskSequence, expansionFilter);
motionMaskSequenceRaw(motionMaskSequenceRaw<motionThreshold )=0;
motionMaskSequenceRaw(motionMaskSequenceRaw>motionThreshold )=1;
motionMaskSequenceRaw = imdilate(motionMaskSequenceRaw, dilateStructure);

% Refine some frames that cannot be correctly detected
motionMaskSequenceRefine = imfilter(motionMaskSequence, expansionFilter);
motionMaskSequenceRefine(motionMaskSequenceRefine<motionThresholdRefine )=0;
motionMaskSequenceRefine(motionMaskSequenceRefine>motionThresholdRefine )=1;
motionMaskSequenceRefine = imdilate(motionMaskSequenceRefine, dilateStructure);

% Replace incorrect frames with refined frames
motionMaskSequence = motionMaskSequenceRaw;
motionMaskSequence(:,:,30:38) =  motionMaskSequenceRefine(:,:,30:38);

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
       
        previousFrameDiff = abs(currentFrame - previousFrame);
        previousFrameDiff2 = abs(currentFrame - previousFrame2);

        blotchMask = zeros(size(currentFrame));
        motionMask = motionMaskSequence(:,:,k);
        
        motionMaskBin = imbinarize(motionMask);
        motionMask = imfill(motionMaskBin,'holes');  
        
        motionPercentage = sum(motionMask(:) == 1)/pixelCount;
        
        if (motionPercentage > 0.75)
           motionMask(:) = 1; 
        end
        
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
        
        %blotchMask = imdilate(blotchMask, strel('disk', 2, 4));
        
        blotchMaskSequence(:,:,k)=blotchMask;

    end 
        
    for i = sceneClipPoints(p)+1:sceneClipPoints(p+1)        
        k = i-first+1;

        %if (k < sceneCutFrames(verticalArtifactSequence))           
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

%              imageSequence(:,:,k) = restoredFrame;
%              imageSequence(:,:,k-1) = restoredPreviousFrame;
            % Recover the intensity gain/loss from averaging pixel values
            imageSequence(:,:,k) = histeq(restoredFrame, imhist(imageSequence(:,:,k)));
            imageSequence(:,:,k-1) = histeq(restoredPreviousFrame, imhist(imageSequence(:,:,k-1)));

        %end
    end   
end

toc

%% Vertical Artifacts Reduction
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

            disp("Done");
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
% for p = 1 : 2 : sceneClipPointsCount-1
%     for i = sceneClipPoints(p)+1:sceneClipPoints(p+1)    
%         k = i-first+1;    
%         
%         currentFrame = uint8(imageSequence(:,:,k));
%         previousFrame = uint8(imageSequence(:,:,k-1));
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
%     end
% end
%%
% Overlay the Text
for f = 1 : length(sceneCutFrames)
    detectedFrame = imageSequence(:,:,sceneCutFrames(f));
    detectedFrame = insertText(detectedFrame,[0,0],'Scene Cut Detected','FontSize',24); 
    imageSequence(:,:,sceneCutFrames(f)) = rgb2gray(detectedFrame);
end

% Save the result
save_sequence(imageSequence, outputPath, prefix, first, digits);

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