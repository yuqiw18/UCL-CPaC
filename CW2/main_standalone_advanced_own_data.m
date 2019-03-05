% If you ran any task related to "gjbLookAtTarget_" image set before please
% type 'clear' in command window to clear the cache before running this
% part! The cache only works for "myData_" image set! (can not be reused with any other section)
clearvars -except flowsFile flowsData imageSequence distanceMatrixAdvanced;
clc;
warning('off','all');

%% Setup
path = 'mydata';
prefix = 'mydata_';
first = 0;
last = 24;
digits = 4;
suffix = 'png';
outputPath = 'output';

%% Initialisation
% Load Optical Flow
if (exist('flowsFile','var') == 0)
    disp("@Load Optical Flow Data");
    tic  
        flowsFile = load('flow/myflows.mat');
        flowsData = flowsFile.flows_a;
    toc
else
    disp("@Optical Flow Already Loaded");
end

%Load Image Sequence
if (exist('imageSequence','var') == 0)
    disp("@Load Image Sequence");
    tic  
        imageSequence = load_sequence_color(path, prefix, first, last, digits, suffix);
    toc
else
    disp("@Image Sequence Already Loaded");
end

imageCount = last - first + 1;

%% Computation
% 1. Compute a Distance Matrix that encodes the similarity in appearance
% between all the frames in the collection.
if (exist('distanceMatrixAdvanced','var') == 0)
    distanceMatrixAdvanced = ComputeDistanceMatrixAdvanced(imageSequence, flowsData);
else
    disp("@Advanced Distance Matrix Already Calculated");
end

[w,h] = size(distanceMatrixAdvanced);

% 2. Convert the Distance Matrix into a graph.
connectionMatrix = DistanceMatrixRejection(distanceMatrixAdvanced, 3);
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

% Define a curve
imshow(imageSequence(:,:,:,selectedImageIndex)),title('Draw a path with at least 5 points');
pointCount = 0;
while(pointCount<5)
    [pathX, pathY]=getline();
    pointCount=size(pathX,1);
    if(pointCount<5)
        % Ask user to draw again if not enough points are collected
        imshow(imageSequence(:,:,:,selectedImageIndex)),title('At least 5 points, please draw again');
    end
end

% Points for reproducing the result in the report
% pointCount = 8; %@frame 2
% pathX = [250;288;304;294;251;238;284;270];
% pathY = [140;145;127;96;89;118;124;107];

% 3. Compute the shortest path for start point.
% In the iteration compute the shortest path for each node in the graph.
[~,paths,~] = graphshortestpath(sparseDistanceMatrix,selectedImageIndex);

% 4,5. Compute the advected location of point s in every other node, using
% optical flow. Pick the path in Paths whose advected location comes closest to the selected point. 

EstimatedClosestAdvLoc = zeros(pointCount,2);
EstimatedClosestAdvLoc(1,:)=[pathX(1),pathY(1)];

disp("@Compute Closest Advected Path");
tic  
    for i = 1:pointCount-1    
        
        startPoint = [pathX(i),pathY(i)];   
        endPoint = [pathX(i+1),pathY(i+1)]; 
        % Calculate the advected position value
        advectedPositionValue = ComputeAdvectedPosition(startPoint, paths, flowsData);
        % Find the closest point
        [closestIndex, closestPoint] = FindClosestPoint(endPoint, advectedPositionValue);
        % Add the point to the estimation list
        EstimatedClosestAdvLoc(i+1,:)= closestPoint;
        % Get the closest path
        closestPaths{i} = paths{closestIndex};
        % Update the path for next iteration
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

% Render this path as the output image sequence for this user-drawn segment.
[outputIndex, outputImageSequence] = ConvertPathsToImageSequence(closestPaths, imageSequence);

% Save the result
outputImageSequence = imresize(outputImageSequence, [300 400]);
save_sequence_color(outputImageSequence,outputPath,'output_adv_od_',0,4);


