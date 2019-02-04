clear;
clc;

path = 'footage';
prefix = 'footage_';
first = 1; % 497
last = 657; % 657
digits = 3;
suffix = 'png';

outputPath = 'output';

% Image Initialisation
rawImageSequence = load_sequence(path, prefix, first, last, digits, suffix);
[vertical,horizontal,sequenceLength]=size(rawImageSequence);
pixelCount = vertical*horizontal;

imageSequence = im2double(rawImageSequence);
imageSequence2 = imageSequence;

% Task 1 variables
sequenceFrameDiff = zeros(size(imageSequence));
sceneCutThreshold = 0.25;
sceneCutFrames = [];
sceneClipPoints = [];

% Task 2 variables
sceneFrameReferenceNumber = [5,7,1];


% Task 3 variables
motionMaskSequence = zeros(size(imageSequence));
blotchMaskSequence = zeros(size(imageSequence));

sceneMotionThreshold = [0.24,0.36,0.6];
sceneBlotchThreshold = [0.08,0.1,0.1];

expansionFilter = fspecial('average', 40);
dilateStructure = strel('disk', 28, 4);

% Task 4 variables
verticalArtifactSequence = 2;

imageSequenceEdge = zeros(size(imageSequence));
sharpenFilter = [0,0,0; 0,2,0; 0,0,0];
meanFilter = [1/9,1/9,1/9;1/9,1/9,1/9;1/9,1/9,1/9];
laplacianFilter = fspecial('laplacian',0);

% Task 5 variables
shakeThreshold = 0.1;

%% Task 1: Scene Cut Detection
disp("@Scene Cut Detection");
tic
for i = 2 : sequenceLength-1

        previousFrameDiff = (sum(abs(imageSequence(:,:,i)-imageSequence(:,:,i-1)),'all'))/pixelCount;
        nextFrameDiff = (sum(abs(imageSequence(:,:,i)-imageSequence(:,:,i+1)),'all'))/pixelCount;          
        transitionDiff = abs(previousFrameDiff - nextFrameDiff);
              
        if (transitionDiff > sceneCutThreshold)     
            sceneCutFrames = [sceneCutFrames, i+first-1];
            disp(strcat(' - Scene Cut Detected @', int2str(i+first-1)));
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

%% Task 2: Global Flicker Reduction
disp("@Global Flicker Reduction");
tic
for p = 1 : 2 : sceneClipPointsCount-1
    s = int8(p/2);   
    for i = sceneClipPoints(p):sceneClipPoints(p+1)    
        f = i-first+1;    
        
        restoredFrame = imageSequence(:,:,f);        
        [frameStart, frameEnd] = FindFrameRange(f, p, sceneClipPoints,sceneFrameReferenceNumber(s));       
        
        averageIntensity = mean(imageSequence(:,:,frameStart:frameEnd),3);
        
        targetHistogram = imhist(averageIntensity);        
        restoredFrame = histeq(restoredFrame, targetHistogram);
        
        imageSequence(:,:,f) = restoredFrame;
    end
end
toc

%% Task 3: Blotch Correction

% disp("@Blotch Correction");
% tic
% 
% motionThreshold = 0.28;
% motionThresholdRefine = 0.9;
% blotchThreshold = 0.1;
% 
% % Motion Mask Generation
% for p = 1 : 2 : sceneClipPointsCount-1  
%     for i = sceneClipPoints(p)+1:sceneClipPoints(p+1)
%     f = i-first+1;
%         sequenceFrameDiff(:,:,f) = abs(imageSequence(:,:,f) - imageSequence(:,:,f-1));   
%     end
%     
%     for i = sceneClipPoints(p):sceneClipPoints(p+1) 
%     f = i-first+1;
%     [frameStart, frameEnd] = FindFrameRange(f, p, sceneClipPoints,2);
%     motionMaskSequence(:,:,f) = sum(sequenceFrameDiff(:,:,frameStart:frameEnd),3);  
%     end  
% end
% 
%     motionMaskSequenceRaw = imfilter(motionMaskSequence, expansionFilter);
% 
% for p = 1 : 2 : sceneClipPointsCount-1
%     s = int8(p/2);
%     frameStart = sceneClipPoints(p);
%     frameEnd = sceneClipPoints(p+1);
% 
%     currentSceneSequenceRaw = motionMaskSequenceRaw(:,:,frameStart:frameEnd);
%     currentSceneSequenceRaw(currentSceneSequenceRaw<sceneMotionThreshold(s))=0;
%     currentSceneSequenceRaw(currentSceneSequenceRaw>sceneMotionThreshold(s))=1;
%     %currentSceneSequenceRaw = imdilate(currentSceneSequenceRaw, dilateStructure);
%     motionMaskSequenceRaw(:,:,frameStart:frameEnd) = currentSceneSequenceRaw;
%     
%     %motionMaskSequenceRefine = motionMaskSequenceRaw(:,:,frameStart:frameEnd);   
% end
% 
% % Refine some frames that cannot be correctly detected
% motionMaskSequenceRefine = motionMaskSequence(:,:,32:36);
% motionMaskSequenceRefine = imfilter(motionMaskSequenceRefine, expansionFilter);
% motionMaskSequenceRefine(motionMaskSequenceRefine<motionThresholdRefine )=0;
% motionMaskSequenceRefine(motionMaskSequenceRefine>motionThresholdRefine )=1;
% %motionMaskSequenceRefine = imdilate(motionMaskSequenceRefine, dilateStructure);
% 
% % Replace incorrect frames with refined frames
% motionMaskSequence = motionMaskSequenceRaw;
% motionMaskSequence(:,:,32:36) =  motionMaskSequenceRefine;
% 
% % motionMaskSequenceBin = imbinarize(motionMaskSequence);
% % motionMaskSequence = imfill(motionMaskSequenceBin,'holes'); 
% 
% % Blotch Mask Generation
% for p = 1 : 2 : sceneClipPointsCount-1
%     s = int8(p/2);
%     for i = sceneClipPoints(p):sceneClipPoints(p+1)       
%         f = i-first+1;
%         
%         currentFrame = imageSequence(:,:,f);
%        
%         if (f - 2 < sceneClipPoints(p))
%             previousFrame = currentFrame;
%             previousFrame2 = currentFrame;
%         else
%             if (f - 1 < sceneClipPoints(p))
%                 previousFrame = currentFrame;
%                 previousFrame2 = imageSequence(:,:,f-1);
%             else
%                 previousFrame = imageSequence(:,:,f-1);
%                 previousFrame2 = imageSequence(:,:,f-2);
%             end
%         end
%        
%         previousFrameDiff = abs(currentFrame - previousFrame);
%         previousFrameDiff2 = abs(currentFrame - previousFrame2);
% 
%         for v = 1: vertical
%             for h = 1: horizontal  
%                 % Mark current pixel as blotch
%                 if (previousFrameDiff(v,h)>blotchThreshold && previousFrameDiff2(v,h)>blotchThreshold)
%                     %blotchMask(v,h) = 1;
%                     blotchMaskSequence(v,h,f) = 1;
%                 end
%                 
%                 % Remove false blotch by checking motion mask
%                 if (motionMaskSequence(v,h,f) == 1)
%                     blotchMaskSequence(v,h,f) = 0;
%                 end
%             end
%         end
%            
% 
%     end   
% end
% 
% %blotchMaskSequence = imdilate(blotchMaskSequence, strel('square', 4));
% 
% for p = 1 : 2 : sceneClipPointsCount-1
%     s = int8(p/2);
%     for i = sceneClipPoints(p)+1:sceneClipPoints(p+1)        
%         f = i-first+1;
% 
%         if (f < sceneCutFrames(verticalArtifactSequence))           
%             [frameStart, frameEnd] = FindFrameRange(f, p, sceneClipPoints,3);
%         
%             % Exclude current frame for blotch correction
%             validFrameCount = frameEnd-frameStart;
% 
%             currentFrame = imageSequence(:,:,f);  
%             restoredFrame = imageSequence(:,:,f);
%             restoredPreviousFrame = imageSequence(:,:,f-1);
%             
%             for v = 1: vertical
%                 for h = 1: horizontal  
%                     if ( blotchMaskSequence(v,h,f) == 1)
%                        % Correct blotch by getting the average value of current
%                        % pixel position from previous and next few frames
%                        % excluding current frame itself (blotch value = error)
%                        restoredFrame(v,h) = (sum(imageSequence(v,h,frameStart:frameEnd),3) - currentFrame(v,h))/(validFrameCount);
% 
%                        % Also correct previous frame because the detection
%                        % compares current frames with previous frames to find
%                        % blotches and previous frame may already have blotches
%                        restoredPreviousFrame(v,h) = (sum(imageSequence(v,h,frameStart:frameEnd),3) - currentFrame(v,h))/(validFrameCount);
%                     end
%                 end
%             end
% 
%              imageSequence(:,:,f) = restoredFrame;
%              imageSequence(:,:,f-1) = restoredPreviousFrame;
%              
%         end
%     end 
% end
% toc

%% Task 4: Vertical Artifacts Reduction
disp("@Vertical Artifacts Reduction");
tic
for p = 1 : 2 : sceneClipPointsCount-1
    for i = sceneClipPoints(p):sceneClipPoints(p+1)    
        
        f = i-first+1;  
        
        currentFrame = imageSequence(:,:,f);       
        currentFrameFrequency = zeros(1,horizontal);
        
        if (f >= sceneCutFrames(verticalArtifactSequence))
            
            for h = 1:horizontal
                for v = 1:vertical
                    currentFrameFrequency(1,h) = currentFrameFrequency(1,h) + currentFrame(v,h);
                end
            end
            
            currentFrameFrequency = currentFrameFrequency/vertical;
            smoothFrameFrequency = medfilt1(currentFrameFrequency,12);
                      
            noiseFrequency = currentFrameFrequency - smoothFrameFrequency;
            
            for h = 1:horizontal
                for v = 1:vertical
                    currentFrame(v,h) = currentFrame(v,h) - noiseFrequency(1,h);
                end
            end
            
            imageSequence(:,:,f) = currentFrame;
                    
%             for v = 1:vertical
%                 for m = 5:-1:1
%                 imageSequence(v,:,f) = medfilt1(imageSequence(v,:,f),m);
%                 end
%             end
% 
%             %Recover features from blurring
%             %Laplacian
%             imageSequenceEdge(:,:,f) = imfilter(imageSequence(:,:,f), laplacianFilter, 'replicate');
%             imageSequence(:,:,f) = imageSequence(:,:,f) - imageSequenceEdge(:,:,f);
% 
%             % Sharpening
%             imageSequence(:,:,k) = imsharpen(imageSequence(:,:,k));
%             imageSequence(:,:,k)= imfilter(imageSequence(:,:,k),sharpenFilter) - imfilter(imageSequence(:,:,k), meanFilter);
            
        end
    end
end
toc

%% Task 5: Camera Shake Calibration
% disp("@Camera Shake Calibration");
% tic
% for p = 1 : 2 : sceneClipPointsCount-1
%     
% %     referenceFrame = imageSequence(:,:,sceneClipPoints(p));
% %     referencePoints = detectFASTFeatures(referenceFrame, 'MinContrast', shakeThreshold);
% %     [referenceFeatures, referencePoints] = extractFeatures(referenceFrame, referencePoints);
%     
%     for i = sceneClipPoints(p)+1:sceneClipPoints(p+1)    
%         k = i-first+1;    
%         
%         currentFrame = imageSequence(:,:,k);
%         previousFrame = imageSequence(:,:,k-1);
%         
%         % detect corners of prev and curr frame
%         currentPoints = detectFASTFeatures(currentFrame, 'MinContrast', shakeThreshold);
%         previousPoints = detectFASTFeatures(previousFrame, 'MinContrast', shakeThreshold);
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
%         tform = estimateGeometricTransform(previousPoints, currentPoints, 'affine');
%         imageSequence(:,:,k) = imwarp(previousFrame, tform, 'OutputView',imref2d(size(currentFrame)));
%         
%     end
% end
% toc

%%
% Overlay the Text
% for t = 1 : length(sceneCutFrames)
%     detectedFrame = imageSequence(:,:,sceneCutFrames(t));
%     detectedFrame = insertText(detectedFrame,[0,0],'Scene Cut Detected','FontSize',24); 
%     imageSequence(:,:,sceneCutFrames(t)) = rgb2gray(detectedFrame);
% end

% Save the result
save_sequence(imageSequence, outputPath, prefix, first, digits);

%save_sequence(motionMaskSequence, 'output2', prefix, first, digits);

% Frame by frame comparasion
%implay([im2double(rawImageSequence), imageSequence]);
implay([blotchMaskSequence,motionMaskSequence,sequenceFrameDiff, im2double(rawImageSequence), imageSequence]);
%implay([im2double(rawImageSequence), imageSequence, abs(im2double(rawImageSequence)-imageSequence)]);
%%
% Functions
% Find closest start and end frames
function [frameStart, frameEnd] = FindFrameRange(currentFrame, currentSequenceNumber, sceneClipPoints, n)
    
startLimit = sceneClipPoints(currentSequenceNumber);
    endLimit = sceneClipPoints(currentSequenceNumber+1);
    
    % Find start frame
    for i = n:-1:0
        if (currentFrame - i < startLimit)
            
        else
            frameStart = currentFrame - i;
            break;
        end
    end
    
    % Find end frame
    for i = n:-1:0
        if (currentFrame + i > endLimit)
            
        else
            frameEnd = currentFrame + i;
            break;
        end
    end
    
end

function [frameStart, frameEnd] = FindFrameRange2(currentFrame, currentSequenceNumber, sceneClipPoints, n)
    
    startLimit = sceneClipPoints(currentSequenceNumber);
    endLimit = sceneClipPoints(currentSequenceNumber+1);
    
    startLimitShift = 0;
    endLimitShift = 0;
    
    % Find start frame
    for i = n:-1:0
        if (currentFrame - i < startLimit)
            
        else
            frameStart = currentFrame - i;
            endLimitShift = n - i;
            break;
        end
    end
    
    % Find end frame
    for i = n:-1:0
        if (currentFrame + i > endLimit)
            
        else
            frameEnd = currentFrame + i;
            startLimitShift = n - i;
            break;
        end
    end
    
    frameStart = frameStart - startLimitShift;
    frameEnd = frameEnd + endLimitShift;
    
end