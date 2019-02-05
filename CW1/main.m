clear;
clc;
warning('off','all');

path = 'footage';
prefix = 'footage_';
first = 1;
last = 657; % 657
digits = 3;
suffix = 'png';

outputPath = 'output';

% Image Initialisation
rawImageSequence = load_sequence(path, prefix, first, last, digits, suffix);
[vertical,horizontal,sequenceLength]=size(rawImageSequence);
pixelCount = vertical*horizontal;

imageSequence = im2double(rawImageSequence);

% % Task 1 variables
% sequenceFrameDiff = zeros(size(imageSequence));
% sceneCutThreshold = 0.25;
% sceneCutFrames = [];
% sceneClipPoints = [];

% % Task 3 variables
% motionMaskSequence = zeros(size(imageSequence));
% blotchMaskSequence = zeros(size(imageSequence));
% 
% sceneMotionThreshold = [0.3,0.3,0.6];
% 
% expansionFilter = fspecial('average', 40);
% dilateStructure = strel('disk', 30, 4);

% % Task 4 variables
% verticalArtifactSequence = 2;
% 
% imageSequenceEdge = zeros(size(imageSequence));
% sharpenFilter = [0,0,0; 0,2,0; 0,0,0];
% meanFilter = [1/9,1/9,1/9;1/9,1/9,1/9;1/9,1/9,1/9];
% laplacianFilter = fspecial('laplacian',0);

%% Task 1: Scene Cut Detection
% disp("@Scene Cut Detection");
% tic
% for i = 2 : sequenceLength-1
% 
%         previousFrameDiff = (sum(abs(imageSequence(:,:,i)-imageSequence(:,:,i-1)),'all'))/pixelCount;
%         nextFrameDiff = (sum(abs(imageSequence(:,:,i)-imageSequence(:,:,i+1)),'all'))/pixelCount;          
%         transitionDiff = abs(previousFrameDiff - nextFrameDiff);
%               
%         if (transitionDiff > sceneCutThreshold)     
%             sceneCutFrames = [sceneCutFrames, i+first-1];
%             disp(strcat(' - Scene Cut Detected @', int2str(i+first-1)));
%         end
% end
% toc
% 
% sceneClipPoints = sceneCutFrames;
% % Append first frame
% sceneClipPoints = [first, sceneClipPoints];
% % Append last frame
% sceneClipPoints = [sceneClipPoints, last];
% % Remove "redundant" cut scene frame
% sceneCutFrames = sceneCutFrames(2:2:length(sceneCutFrames));

[sceneCutFrames, sceneClipPoints] = SceneCutDetection(first, last, imageSequence);
imageSequence = GlobalFlickerReduction(first, last, imageSequence, sceneClipPoints);
imageSequence = BlotchCorrection(first, last, imageSequence, sceneCutFrames, sceneClipPoints);
imageSequence = VerticalArtifactReduction(first, last, imageSequence, sceneCutFrames, sceneClipPoints);
imageSequence = CameraShakeCalibration(first, last, imageSequence, sceneClipPoints);
implay([im2double(rawImageSequence), imageSequence]);

sceneClipPointsCount = length(sceneClipPoints);

%% Task 2: Global Flicker Reduction
% disp("@Global Flicker Reduction");
% tic
% for p = 1 : 2 : sceneClipPointsCount-1
%     s = int8(p/2);   
%     for i = sceneClipPoints(p):sceneClipPoints(p+1)    
%         f = i-first+1;    
%         
%         restoredFrame = imageSequence(:,:,f);        
%         [frameStart, frameEnd] = FindFrameRange(f, p, sceneClipPoints,6);       
%         
%         averageIntensity = mean(imageSequence(:,:,frameStart:frameEnd),3);
%         
%         targetHistogram = imhist(averageIntensity);        
%         restoredFrame = histeq(restoredFrame, targetHistogram);
%         
%         imageSequence(:,:,f) = restoredFrame;
%     end
% end
% toc

%% Task 3: Blotch Correction
% disp("@Blotch Correction");
% tic
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
%     
%     frameStart = sceneClipPoints(p);
%     frameEnd = sceneClipPoints(p+1);
% 
%     currentSceneSequenceRaw = motionMaskSequenceRaw(:,:,frameStart:frameEnd);
%     currentSceneSequenceRaw(currentSceneSequenceRaw<sceneMotionThreshold(s))=0;
%     currentSceneSequenceRaw(currentSceneSequenceRaw>sceneMotionThreshold(s))=1;
%     motionMaskSequenceRaw(:,:,frameStart:frameEnd) = currentSceneSequenceRaw;
%      
% end
% 
% % %motionMaskSequenceRaw = imdilate(motionMaskSequenceRaw, dilateStructure);
% % 
% % % Refine some frames that cannot be correctly detected
% % motionMaskSequenceRefine = motionMaskSequence(:,:,32:36);
% % motionMaskSequenceRefine = imfilter(motionMaskSequenceRefine, expansionFilter);
% % motionMaskSequenceRefine(motionMaskSequenceRefine<motionThresholdRefine )=0;
% % motionMaskSequenceRefine(motionMaskSequenceRefine>motionThresholdRefine )=1;
% % %motionMaskSequenceRefine = imdilate(motionMaskSequenceRefine, dilateStructure);
% % 
% % % Replace incorrect frames with refined frames
% motionMaskSequence = motionMaskSequenceRaw;
% % motionMaskSequence(:,:,32:36) =  motionMaskSequenceRefine;
% 
% motionMaskSequence = imdilate(motionMaskSequence, dilateStructure);
% %motionMaskSequence = imfill(motionMaskSequenceBin,'holes'); 
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
% disp("@Vertical Artifacts Reduction");
% tic
% for p = 1 : 2 : sceneClipPointsCount-1
%     for i = sceneClipPoints(p):sceneClipPoints(p+1)    
%         
%         f = i-first+1;  
%         
%         currentFrame = imageSequence(:,:,f);       
%         currentFrameFrequency = zeros(1,horizontal);
%         
%         if (f >= sceneCutFrames(verticalArtifactSequence))
%             
% %             for h = 1:horizontal
% %                 for v = 1:vertical
% %                     currentFrameFrequency(1,h) = currentFrameFrequency(1,h) + currentFrame(v,h);
% %                 end
% %             end
% %             
% %             currentFrameFrequency = currentFrameFrequency/vertical;
% %             smoothFrameFrequency = medfilt1(currentFrameFrequency,9);
% %                       
% %             noiseFrequency = currentFrameFrequency - smoothFrameFrequency;
% %             
% %             for h = 1:horizontal
% %                 for v = 1:vertical
% %                     currentFrame(v,h) = currentFrame(v,h) - noiseFrequency(1,h);
% %                 end
% %             end
% %             
% %             imageSequence(:,:,f) = currentFrame;
%                     
%             for v = 1:vertical
%                 %for m = 1:1:5
%                 imageSequence(v,:,f) = medfilt1(imageSequence(v,:,f),5);
%                 %end
%             end
% 
%             %Recover features from blurring
%             %Laplacian
%             imageSequenceEdge(:,:,f) = imfilter(imageSequence(:,:,f), laplacianFilter, 'replicate');
%             imageSequence(:,:,f) = imageSequence(:,:,f) - imageSequenceEdge(:,:,f);
% 
% %             % Sharpening
% %             imageSequence(:,:,f) = imsharpen(imageSequence(:,:,f));
% %             imageSequence(:,:,f)= imfilter(imageSequence(:,:,f),sharpenFilter) - imfilter(imageSequence(:,:,f), meanFilter);
% 
% %             smoothFrame = imageSequence(:,:,f);
% %             smoothFrameFrequency = zeros(1,horizontal);
% %             for h = 1:horizontal
% %                 for v = 1:vertical
% %                     smoothFrameFrequency(1,h) = smoothFrameFrequency(1,h) + smoothFrame(v,h);
% %                 end
% %             end
% %             
% %             figure;
% %             subplot(2,1,1);
% %             plot(currentFrameFrequency);
% %             subplot(2,1,2);
% %             plot(smoothFrameFrequency);
%         end
%     end
% end
% toc

%% Task 5: Camera Shake Calibration




% disp("@Camera Shake Calibration");
% tic
% for p = 1 : 2 : sceneClipPointsCount-1
%     
%     % Initialisation
%     referenceFrame = mean(imageSequence(:,:,sceneClipPoints(p):sceneClipPoints(p+1)),3);
%     [row, col, ~] = size(referenceFrame);
%     worldCoord = imref2d(size(referenceFrame));
%     windowRow = 2 * row -1;
%     windowCol = 2 * col -1;
%     shiftRow = (1 + (windowRow-1)/2); 
%     shiftCol = (1 + (windowCol-1)/2);
%     
%     Ga = fft2(referenceFrame,windowRow,windowCol);
%        
%     for i = sceneClipPoints(p):sceneClipPoints(p+1)    
%         f = i-first+1; 
%         
%         currentFrame = single(imageSequence(:,:,f));      
%         
%         % Phase correlation
%         Gb = fft2(currentFrame,windowRow,windowCol);
%         GbConj = conj(Gb);
%         R = (Ga .* GbConj)./abs(Ga .* GbConj); 
%         r = real(ifft2(R));
%         r = fftshift(r);
%         
%         % Find weighted position of the peak 
%         [peakRow, peakCol] = findWeightedPeak(r, 2);
%         
%         % Compute the shift values
%         frameShiftRow = peakRow-shiftRow;
%         frameShiftCol = peakCol-shiftCol;
%         
%         % Construct transform matrix
%         frameShift = affine2d([1, 0, 0; 0, 1, 0; frameShiftCol, frameShiftRow, 1]);
%         
%         restoredFrame = imwarp(currentFrame, frameShift, 'OutputView',worldCoord);        
%         imageSequence(:,:,f) = restoredFrame;
%              
%     end
% end
% toc

%%
imageSequence = OverlayText(imageSequence, sceneCutFrames);

%Overlay the Text
% for t = 1 : length(sceneCutFrames)
%     detectedFrame = imageSequence(:,:,sceneCutFrames(t));
%     detectedFrame = insertText(detectedFrame,[0,0],'Scene Cut Detected','FontSize',24); 
%     imageSequence(:,:,sceneCutFrames(t)) = rgb2gray(detectedFrame);
% end

% Save the result
%save_sequence(imageSequence, outputPath, prefix, first, digits);

% Frame by frame comparasion
%implay([im2double(rawImageSequence), imageSequence]);
%implay([blotchMaskSequence,motionMaskSequence,sequenceFrameDiff, im2double(rawImageSequence), imageSequence]);


%% Task 1: Scene Cut Detection
function [sceneCutFrames, sceneClipPoints] = SceneCutDetection(first, last, imageSequence)
    disp("@Scene Cut Detection");
    tic  
    
    threshold = 0.25;
    sceneCutFrames = [];
    
    [vertical,horizontal,sequenceLength]=size(imageSequence);
    pixelCount = vertical*horizontal;
    
    for i = 2 : sequenceLength-1

        previousFrameDiff = (sum(abs(imageSequence(:,:,i)-imageSequence(:,:,i-1)),'all'))/pixelCount;
        nextFrameDiff = (sum(abs(imageSequence(:,:,i)-imageSequence(:,:,i+1)),'all'))/pixelCount;          
        transitionDiff = abs(previousFrameDiff - nextFrameDiff);

        if (transitionDiff > threshold)     
            sceneCutFrames = [sceneCutFrames, i+first-1];
            disp(strcat(' - Scene Cut Detected @', int2str(i+first-1)));
        end
    end
    
    sceneClipPoints = sceneCutFrames;
    
    % Append first frame
    sceneClipPoints = [first, sceneClipPoints];
    
    % Append last frame
    sceneClipPoints = [sceneClipPoints, last];
    
    % Remove "redundant" cut scene frame
    sceneCutFrames = sceneCutFrames(2:2:length(sceneCutFrames));

    toc
end

%% Task 2: Global Flicker Reduction
function imageSequence = GlobalFlickerReduction(first, last, imageSequence, sceneClipPoints)
    disp("@Global Flicker Reduction");
    tic
    sceneClipPointsCount = length(sceneClipPoints);
    for p = 1 : 2 : sceneClipPointsCount-1
        for i = sceneClipPoints(p):sceneClipPoints(p+1)    
            f = i-first+1;    

            restoredFrame = imageSequence(:,:,f);        
            [frameStart, frameEnd] = FindFrameRange(f, p, sceneClipPoints,6);       

            averageIntensity = mean(imageSequence(:,:,frameStart:frameEnd),3);

            targetHistogram = imhist(averageIntensity);        
            restoredFrame = histeq(restoredFrame, targetHistogram);

            imageSequence(:,:,f) = restoredFrame;
        end
    end
    toc
end

%% Task 3: Blotch Correction
function imageSequence = BlotchCorrection(first, last, imageSequence, sceneCutFrames, sceneClipPoints)
    disp("@Blotch Correction");
    tic
    [vertical,horizontal,~]=size(imageSequence);
    sceneClipPointsCount = length(sceneClipPoints);
    
    sequenceFrameDiff = zeros(size(imageSequence));
    motionMaskSequence = zeros(size(imageSequence));
    blotchMaskSequence = zeros(size(imageSequence));

    sceneMotionThreshold = [0.3,0.3,0.6];
    blotchThreshold = 0.1;
    
    expansionFilter = fspecial('average', 40);
    dilateStructure = strel('disk', 30, 4);
 
    % Motion Mask Generation
    for p = 1 : 2 : sceneClipPointsCount-1  
        for i = sceneClipPoints(p)+1:sceneClipPoints(p+1)
        f = i-first+1;
            sequenceFrameDiff(:,:,f) = abs(imageSequence(:,:,f) - imageSequence(:,:,f-1));   
        end

        for i = sceneClipPoints(p):sceneClipPoints(p+1) 
        f = i-first+1;
        [frameStart, frameEnd] = FindFrameRange(f, p, sceneClipPoints,2);
        motionMaskSequence(:,:,f) = sum(sequenceFrameDiff(:,:,frameStart:frameEnd),3);  
        end  
    end

        motionMaskSequenceRaw = imfilter(motionMaskSequence, expansionFilter);

    for p = 1 : 2 : sceneClipPointsCount-1
        s = int8(p/2);

        frameStart = sceneClipPoints(p);
        frameEnd = sceneClipPoints(p+1);

        currentSceneSequenceRaw = motionMaskSequenceRaw(:,:,frameStart:frameEnd);
        currentSceneSequenceRaw(currentSceneSequenceRaw<sceneMotionThreshold(s))=0;
        currentSceneSequenceRaw(currentSceneSequenceRaw>sceneMotionThreshold(s))=1;
        motionMaskSequenceRaw(:,:,frameStart:frameEnd) = currentSceneSequenceRaw;

    end

    % %motionMaskSequenceRaw = imdilate(motionMaskSequenceRaw, dilateStructure);
    % 
    % % Refine some frames that cannot be correctly detected
    % motionMaskSequenceRefine = motionMaskSequence(:,:,32:36);
    % motionMaskSequenceRefine = imfilter(motionMaskSequenceRefine, expansionFilter);
    % motionMaskSequenceRefine(motionMaskSequenceRefine<motionThresholdRefine )=0;
    % motionMaskSequenceRefine(motionMaskSequenceRefine>motionThresholdRefine )=1;
    % %motionMaskSequenceRefine = imdilate(motionMaskSequenceRefine, dilateStructure);
    % 
    % % Replace incorrect frames with refined frames
    motionMaskSequence = motionMaskSequenceRaw;
    % motionMaskSequence(:,:,32:36) =  motionMaskSequenceRefine;

    motionMaskSequence = imdilate(motionMaskSequence, dilateStructure);
    %motionMaskSequence = imfill(motionMaskSequenceBin,'holes'); 

    % Blotch Mask Generation
    for p = 1 : 2 : sceneClipPointsCount-1
        s = int8(p/2);
        for i = sceneClipPoints(p):sceneClipPoints(p+1)       
            f = i-first+1;

            currentFrame = imageSequence(:,:,f);

            if (f - 2 < sceneClipPoints(p))
                previousFrame = currentFrame;
                previousFrame2 = currentFrame;
            else
                if (f - 1 < sceneClipPoints(p))
                    previousFrame = currentFrame;
                    previousFrame2 = imageSequence(:,:,f-1);
                else
                    previousFrame = imageSequence(:,:,f-1);
                    previousFrame2 = imageSequence(:,:,f-2);
                end
            end

            previousFrameDiff = abs(currentFrame - previousFrame);
            previousFrameDiff2 = abs(currentFrame - previousFrame2);

            for v = 1: vertical
                for h = 1: horizontal  
                    % Mark current pixel as blotch
                    if (previousFrameDiff(v,h)>blotchThreshold && previousFrameDiff2(v,h)>blotchThreshold)
                        %blotchMask(v,h) = 1;
                        blotchMaskSequence(v,h,f) = 1;
                    end

                    % Remove false blotch by checking motion mask
                    if (motionMaskSequence(v,h,f) == 1)
                        blotchMaskSequence(v,h,f) = 0;
                    end
                end
            end


        end   
    end

    %blotchMaskSequence = imdilate(blotchMaskSequence, strel('square', 4));

    for p = 1 : 2 : sceneClipPointsCount-1
        s = int8(p/2);
        for i = sceneClipPoints(p)+1:sceneClipPoints(p+1)        
            f = i-first+1;

            if (f < sceneCutFrames(2))           
                [frameStart, frameEnd] = FindFrameRange(f, p, sceneClipPoints,3);

                % Exclude current frame for blotch correction
                validFrameCount = frameEnd-frameStart;

                currentFrame = imageSequence(:,:,f);  
                restoredFrame = imageSequence(:,:,f);
                restoredPreviousFrame = imageSequence(:,:,f-1);

                for v = 1: vertical
                    for h = 1: horizontal  
                        if ( blotchMaskSequence(v,h,f) == 1)
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

                 imageSequence(:,:,f) = restoredFrame;
                 imageSequence(:,:,f-1) = restoredPreviousFrame;

            end
        end 
    end
    toc


end

%% Task 4: Vertical Artifact Reduction
function imageSequence = VerticalArtifactReduction(first, last, imageSequence, sceneCutFrames, sceneClipPoints)
    disp("@Vertical Artifacts Reduction");
    tic
    
    imageSequenceEdge = zeros(size(imageSequence));
    sharpenFilter = [0,0,0; 0,2,0; 0,0,0];
    meanFilter = [1/9,1/9,1/9;1/9,1/9,1/9;1/9,1/9,1/9];
    laplacianFilter = fspecial('laplacian',0);
 
    [vertical,horizontal,~]=size(imageSequence);
    sceneClipPointsCount = length(sceneClipPoints);
    for p = 1 : 2 : sceneClipPointsCount-1
        for i = sceneClipPoints(p):sceneClipPoints(p+1)    

            f = i-first+1;  

            currentFrame = imageSequence(:,:,f);       
            currentFrameFrequency = zeros(1,horizontal);

            if (f >= sceneCutFrames(2))

    %             for h = 1:horizontal
    %                 for v = 1:vertical
    %                     currentFrameFrequency(1,h) = currentFrameFrequency(1,h) + currentFrame(v,h);
    %                 end
    %             end
    %             
    %             currentFrameFrequency = currentFrameFrequency/vertical;
    %             smoothFrameFrequency = medfilt1(currentFrameFrequency,9);
    %                       
    %             noiseFrequency = currentFrameFrequency - smoothFrameFrequency;
    %             
    %             for h = 1:horizontal
    %                 for v = 1:vertical
    %                     currentFrame(v,h) = currentFrame(v,h) - noiseFrequency(1,h);
    %                 end
    %             end
    %             
    %             imageSequence(:,:,f) = currentFrame;

                for v = 1:vertical
                    %for m = 1:1:5
                    imageSequence(v,:,f) = medfilt1(imageSequence(v,:,f),5);
                    %end
                end

                %Recover features from blurring
                %Laplacian
                imageSequenceEdge(:,:,f) = imfilter(imageSequence(:,:,f), laplacianFilter, 'replicate');
                imageSequence(:,:,f) = imageSequence(:,:,f) - imageSequenceEdge(:,:,f);

    %             % Sharpening
    %             imageSequence(:,:,f) = imsharpen(imageSequence(:,:,f));
    %             imageSequence(:,:,f)= imfilter(imageSequence(:,:,f),sharpenFilter) - imfilter(imageSequence(:,:,f), meanFilter);

    %             smoothFrame = imageSequence(:,:,f);
    %             smoothFrameFrequency = zeros(1,horizontal);
    %             for h = 1:horizontal
    %                 for v = 1:vertical
    %                     smoothFrameFrequency(1,h) = smoothFrameFrequency(1,h) + smoothFrame(v,h);
    %                 end
    %             end
    %             
    %             figure;
    %             subplot(2,1,1);
    %             plot(currentFrameFrequency);
    %             subplot(2,1,2);
    %             plot(smoothFrameFrequency);
            end
        end
    end
    toc
end


%% Task 5: Camera Shake Calibration
function imageSequence = CameraShakeCalibration(first, last, imageSequence, sceneClipPoints)
    disp("@Camera Shake Calibration");
    tic
    sceneClipPointsCount = length(sceneClipPoints);
    for p = 1 : 2 : sceneClipPointsCount-1

        % Initialisation
        referenceFrame = mean(imageSequence(:,:,sceneClipPoints(p):sceneClipPoints(p+1)),3);
        [row, col, ~] = size(referenceFrame);
        worldCoord = imref2d(size(referenceFrame));
        windowRow = 2 * row -1;
        windowCol = 2 * col -1;
        shiftRow = (1 + (windowRow-1)/2); 
        shiftCol = (1 + (windowCol-1)/2);

        Ga = fft2(referenceFrame,windowRow,windowCol);

        for i = sceneClipPoints(p):sceneClipPoints(p+1)    
            f = i-first+1; 

            currentFrame = single(imageSequence(:,:,f));      

            % Phase correlation
            Gb = fft2(currentFrame,windowRow,windowCol);
            GbConj = conj(Gb);
            R = (Ga .* GbConj)./abs(Ga .* GbConj); 
            r = real(ifft2(R));
            r = fftshift(r);

            % Find weighted position of the peak 
            [peakRow, peakCol] = findWeightedPeak(r, 2);

            % Compute the shift values
            frameShiftRow = peakRow-shiftRow;
            frameShiftCol = peakCol-shiftCol;

            % Construct transform matrix
            frameShift = affine2d([1, 0, 0; 0, 1, 0; frameShiftCol, frameShiftRow, 1]);

            restoredFrame = imwarp(currentFrame, frameShift, 'OutputView',worldCoord);        
            imageSequence(:,:,f) = restoredFrame;

        end
    end
    toc
end

% Compute weighted peak position: subpixel precision
function [weightedPeakRow, weightedPeakCol] = findWeightedPeak(matrix, span)

    % Intialise
    [rowLimit, colLimit] = size(matrix);
    weightedSum = 0;
    weightedCol = 0;
    weightedRow = 0;

    % Find the peak value and its position
    peakValue = max(matrix(:));
    [peakRow,peakCol] = find(matrix==peakValue);

    % Setup the (1+span)x(1+span) area
    startRow = peakRow - span;
    endRow = peakRow + span;
    startCol = peakCol - span;
    endCol = peakCol + span;

    % Clipping if out of boundaries
    if (startRow < 0)
        startRow = 0;
    end
    if (endRow > rowLimit)
        endRow = rowLimit;
    end
    if (startCol < 0)
        startCol = 0;
    end
    if (endCol > colLimit)
        endCol = colLimit;
    end

    % Compute the weighted sum
    for r=startRow:endRow
       for c = startCol:endCol
            weightedCol = weightedCol + c * matrix(r,c);
            weightedRow = weightedRow + r * matrix(r,c);
            weightedSum = weightedSum + matrix(r,c);
       end 
    end

    % Assign values
    weightedPeakCol = weightedCol/weightedSum;
    weightedPeakRow = weightedRow/weightedSum;

end

function imageSequence = OverlayText(imageSequence, sceneCutFrames)
    for t = 1 : length(sceneCutFrames)
        detectedFrame = imageSequence(:,:,sceneCutFrames(t));
        detectedFrame = insertText(detectedFrame,[0,0],'Scene Cut Detected','FontSize',24); 
        imageSequence(:,:,sceneCutFrames(t)) = rgb2gray(detectedFrame);
    end
end

% Find closest start and end frames
function [frameStart, frameEnd] = FindFrameRange(currentFrame, currentSequenceNumber, sceneClipPoints, n)
    
startLimit = sceneClipPoints(currentSequenceNumber);
    endLimit = sceneClipPoints(currentSequenceNumber+1);
    
    % Find the start frame
    for i = n:-1:0
        if (currentFrame - i < startLimit)            
        else
            frameStart = currentFrame - i;
            break;
        end
    end
    
    % Find the end frame
    for i = n:-1:0
        if (currentFrame + i > endLimit)          
        else
            frameEnd = currentFrame + i;
            break;
        end
    end 
end
