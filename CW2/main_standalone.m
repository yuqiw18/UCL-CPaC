% If you ran any task related to "myData_" image set before please
% type 'clear' in command window to clear the cache before running this
% part! The cache only works for "gjbLookAtTarget_" image set! (can only be reused with its advanced section)
clearvars -except flowsFile flowsData imageSequence distanceMatrix distanceMatrixAdvanced;
clc;
warning('off','all');

% Setup
path = 'data';
prefix = 'gjbLookAtTarget_';
first = 0;
last = 71;
digits = 4;
suffix = 'jpg';
outputPath = 'output';

% Initialisation
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
if (exist('distanceMatrix','var') == 0)
    distanceMatrix = ComputeDistanceMatrix(imageSequence);
else
    disp("@Distance Matrix Already Calculated");
end

% 2. Convert the Distance Matrix into a graph.
connectionMatrix = DistanceMatrixRejection(distanceMatrix, 7);
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
% pointCount = 7; %@frame 3
% pathX = [243,244,229,259,229,254,253];
% pathY = [191,176,171,170,166,158,192];

% pointCount = 6; %@frame 53
% pathX = [274;256;245;246;263;284];
% pathY = [167;169;179;195;199;195];

% 3. Compute the shortest path for start point.
% In the iteration compute the shortest path for each node in the graph.
[~,paths,~] = graphshortestpath(sparseDistanceMatrix,selectedImageIndex);

% 4,5. Compute the advected location of point s in every other node, using
% optical flow. Pick the path in Paths whose advected location comes closest to the selected point. 

EstimatedClosestAdvLoc = zeros(pointCount-1,2);
EstimatedClosestAdvLoc(1,:)=[pathX(1),pathY(1)];

disp("@Compute Closest Advected Path");
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

% Render this path as the output image sequence for this user-drawn segment.
[~, outputImageSequence] = ConvertPathsToImageSequence(closestPaths, imageSequence);

% Save the result
outputImageSequence = imresize(outputImageSequence, [300 400]);
save_sequence_color(outputImageSequence,outputPath,'output_basic_',0,4);
