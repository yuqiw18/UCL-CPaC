clearvars -except flowsFile flowsData imageSequence distanceMatrixAdvanced distanceMatrix;
clc;
warning('off','all');

%% Setup
path = 'data';
prefix = 'gjbLookAtTarget_';
first = 0;
last = 71;
digits = 4;
suffix = 'jpg';
outputPath = 'output';

%% Initialisation
% Load Optical Flow
if (exist('flowsFile','var') == 0)
    disp("@Loading Optical Flow Data");
    tic  
        flowsFile = load('flow/flows.mat');
        flowsData = flowsFile.flows_a;
    toc
else
    disp("@Optical Flow Already Loaded");
end

%Load Image Sequence
if (exist('imageSequence','var') == 0)
    disp("@Loading Image Sequence");
    tic  
        rawImageSequence = load_sequence_color(path, prefix, first, last, digits, suffix);
        % Downscale image size to 30% to match optical flow used in later stage
        imageSequence = imresize(rawImageSequence, 0.3);
    toc
else
    disp("@Image Sequence Already Loaded");
end

imageCount = last - first + 1;

%% Basic Section
% 1. Compute a Distance Matrix that encodes the similarity in appearance
% between all the frames in the collection.
if (exist('distanceMatrixAdvanced','var') == 0)
    distanceMatrixAdvanced = ComputeDistanceMatrixAdvanced(imageSequence, flowsData);
else
    disp("@Advanced Distance Matrix Already Calculated");
end

[w,h] = size(distanceMatrixAdvanced);

% 2. Convert the Distance Matrix into a graph.
connectionMatrix = DistanceMatrixRejection(distanceMatrixAdvanced);
sparseDistanceMatrix = sparse(connectionMatrix);
% graph = biograph(sparseDistanceMatrix,[],'ShowArrows','off','LayoutType','equilibrium');
% view(graph);

% User input
% Specify the start frame
selectedImageIndex = input('Please specify the frame index to start:');
if (selectedImageIndex < first || selectedImageIndex > last)
    disp('@Invalid Frame: Default - First Frame')
    selectedImageIndex = first+1;
end

imshow(imageSequence(:,:,:,selectedImageIndex)),title('Draw a path with at least 5 points');

% pointCount = 0;
% 
% while(pointCount<5)
%     [pathX, pathY]=getline();
%     pointCount=size(pathX,1);
%     if(pointCount<5)
%         % Ask user to draw again if not enough points are collected
%         imshow(imageSequence(:,:,:,selectedImageIndex)),title('At least 5 points, please draw again');
%     end
% end

pointCount = 6;
pathX = [274.9535073409462,237.3678629690049,224.8393148450244,227.1884176182708,247.5473083197390,289.8311582381729];
pathY = [151.9192495921696,154.2683523654159,169.1460032626427,189.5048939641109,206.7316476345840,204.3825448613376];

% 3. For each node in the graph, compute the shortest path.
[~,paths,~] = graphshortestpath(sparseDistanceMatrix,selectedImageIndex);

% 4. Compute the advected location of point s in every other node, using
% optical flow.

EstimatedClosestAdvLoc = zeros(pointCount,2);
EstimatedClosestAdvLoc(1,:)=[pathX(1),pathY(1)];

disp("@Computing Closest Advected Path");
tic  
    for i = 1:pointCount-1    
        startPoint = [pathX(i),pathY(i)];   
        endPoint = [pathX(i+1),pathY(i+1)]; 
        advectedPositionValue = ComputeAdvectedPosition(startPoint, paths, flowsData);      
        [closestIndex, closestPoint] = FindClosestPoint(endPoint, advectedPositionValue); 
        EstimatedClosestAdvLoc(i+1,:)= closestPoint;   
        closestPaths{i} = paths{closestIndex};   
        [~,paths,~] = graphshortestpath(sparseDistanceMatrix,closestIndex);
    end
toc

% Show selected points and estimated points
close all;
figure;
imshow(imageSequence(:,:,:,selectedImageIndex));
hold on;
plot(pathX,pathY,'g*-');
hold on;
plot(EstimatedClosestAdvLoc(:,1), EstimatedClosestAdvLoc(:,2), 'ro-');
hold off;

% 5. Pick the path in Paths whose advected location comes closest to the selected point. 
% Render this path as the output image sequence for this user-drawn segment.
[outputIndex, outputImageSequence] = ConvertPathsToImageSequence(closestPaths, imageSequence);

%Interpolated frames
frame = 6;
outputIndexInterpolated = frame * (length(outputIndex)-1) + length(outputIndex);
% e.g. F-123456-F-123456-F-123456-F
% F raw frame + frame * (F-1)

[height,width,channel,~] = size(outputImageSequence);
interpolatedImageSequence = zeros(height,width,channel,outputIndexInterpolated);

disp("@Performing Motion Interpolation");
    tic  
    for i = 1: length(outputIndex)-1

        currentFrameIndex = outputIndex(i);
        nextFrameIndex = outputIndex(i+1);

        currentFrame = imageSequence(:,:,:,currentFrameIndex);
        nextFrame = imageSequence(:,:,:,nextFrameIndex);

        % Get the optical flow between the current frame and the next frame
        % Since the previous step has removed duplicated connecting frames
        % currentFrameIndex will never be the same as nextFrameIndex
        if (currentFrameIndex > nextFrameIndex)
            k = (currentFrameIndex-1)*(currentFrameIndex-2)/2+nextFrameIndex;
            flow = flowsData(:,:,:,k);
        else
            k = (nextFrameIndex-1)*(nextFrameIndex-2)/2+currentFrameIndex;
            flow = -flowsData(:,:,:,k);
        end

        interpolatedImage = MotionInterpolation(currentFrame,nextFrame,flow,frame);
        interpolatedImageSequence(:,:,:,(i-1)*(frame+1)+1) = currentFrame;
        for f = 1:frame
            interpolatedImageSequence(:,:,:,(i-1)*(frame+1)+1+f) = interpolatedImage(:,:,:,f);
        end    
    end
    interpolatedImageSequence(:,:,:,outputIndexInterpolated) = imageSequence(:,:,:,outputIndex(length(outputIndex)));
toc

implay(interpolatedImageSequence);
save_sequence_color(interpolatedImageSequence,outputPath,'output_adv_',0,4);


