clearvars -except flowsFile flowsData imageSequence distanceMatrix distanceMatrixAdvanced;
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
    disp("@Load Optical Flow Data");
    tic  
        flowsFile = load('flow/flows.mat');
        flowsData = flowsFile.flows_a;
    toc
else
    disp("@Optical Flow Already Loaded");
end

%Load Image Sequence
if (exist('imageSequence','var') == 0)
    disp("@Load Image Sequence");
    tic  
        rawImageSequence = load_sequence_color(path, prefix, first, last, digits, suffix);
        imageSequence = imresize(rawImageSequence, 0.3);
    toc
else
    disp("@Image Sequence Already Loaded");
end
% Downsize image size to 30% to match optical flow used in later stage

imageCount = last - first + 1;
% imageCount = size(imageSequence,4);

%% Basic Section
% 1. Compute a Distance Matrix that encodes the similarity in appearance
% between all the frames in the collection.
if (exist('distanceMatrix','var') == 0)
    distanceMatrix = ComputeDistanceMatrix(imageSequence);
else
    disp("@Distance Matrix Already Calculated");
end

% 2. Convert the Distance Matrix into a graph.

connectionMatrix = DistanceMatrixRejection(distanceMatrix);
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

% 3. For each node in the graph, compute the shortest path.
[~,paths,~] = graphshortestpath(sparseDistanceMatrix,selectedImageIndex);
[S, ~] = graphconncomp(sparseDistanceMatrix, 'Weak', true);

% 4. Compute the advected location of point s in every other node, using
% optical flow.

EstimatedClosestAdvLoc = zeros(pointCount-1,2);
EstimatedClosestAdvLoc(1,:)=[pathX(1),pathY(1)];

disp("@Compute Closest Advected Path");
tic  
    for i = 1:pointCount-1 
        
%         % Swap X, Y 
%         startPoint = [pathY(i),pathX(i)];   
%         endPoint = [pathY(i+1),pathX(i+1)];    
%         [closestIndex, closestX, closestY] = FindShortestAdvectedPath(startPoint, endPoint, paths, flowsData);
%         EstimatedClosestAdvLoc(i+1,:)= [closestX, closestY];
        
        startPoint = [pathX(i),pathY(i)];   
        endPoint = [pathX(i+1),pathY(i+1)]; 
        
        advectedPositionValue = ComputeAdvectedPosition(startPoint, paths, flowsData); 
        
        [closestIndex, closestPoint] = FindClosestPoint(endPoint, advectedPositionValue); 
        
        EstimatedClosestAdvLoc(i+1,:)= closestPoint;  
        
        closestPaths{i} = paths{closestIndex};   
        
        [~,paths,~] = graphshortestpath(sparseDistanceMatrix,closestIndex);
    end
toc

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
outputImageSequence = ConvertPathsToImageSequence(closestPaths, imageSequence);
implay(outputImageSequence);

save_sequence_color(outputImageSequence,outputPath,'output_',0,4);


