clearvars -except flowsFile flowsData imageSequence distanceMatrix;
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
if exist('flowsFile','var') == 0
    disp("@Load Optical Flow Data");
    tic  
        flowsFile = load('flow/flows.mat');
        flowsData = flowsFile.flows_a;
    toc
else
    disp("@Optical Flow Already Loaded");
end

%Load Image Sequence
if exist('imageSequence','var') == 0
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

if exist('distanceMatrix','var') == 0
    distanceMatrix = ComputeDistanceMatrix(imageSequence);
else
    disp("@Distance Matrix Already Calculated");
end

% figure;
% title("Distance Matrix");
% imshow(distanceMatrix);

% probabilityMatrix = ComputeProbabilityMatrix(distanceMatrix);
% figure;
% title("Probability Matrix");
% imshow(probabilityMatrix);

[w,h] = size(distanceMatrix);
connectionMatrix = distanceMatrix;

% 2. Convert the Distance Matrix into a graph.
for i = 1:w-1   
    average = sum(connectionMatrix(i,i:w),'all')/(w-i);  
    %disp(average);
    for j = i:w
        if (connectionMatrix(i,j)>average)
            connectionMatrix(i,j)=0;
        end
        if (connectionMatrix(j,i)>average)
            connectionMatrix(j,i)=0;
        end
    end
end

sparseDistanceMatrix = sparse(distanceMatrix);
% graph = biograph(sparseDistanceMatrix,[],'ShowArrows','off','LayoutType','equilibrium');
% view(graph);

% graph = graphminspantree(sparseDistanceMatrix);

% User input
selectedImageIndex = input('Please specify the frame index to start:');
if (selectedImageIndex < first || selectedImageIndex > last)
    disp('@Invalid Frame: Default - First Frame')
    selectedImageIndex = first+1;
end

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

% 4. Compute the advected location of point s in every other node, using
% optical flow.

EstimatedClosestAdvLoc = zeros(pointCount,2);
% EstimatedClosestAdvLoc = [pathX, pathY]
EstimatedClosestAdvLoc(1,:)=[pathX(1),pathY(1)];
imageSequenceIndex = [];

disp("@Compute Shortest Path");
tic  
    for i = 1:pointCount
        currentAdvLoc = [pathY(i),pathX(i)];
        %nextAdvLoc = [pathY(i+1),pathX(i+1)]; 
        if (i+1>pointCount)
           nextAdvLoc = [pathY(i),pathX(i)];
        else
           nextAdvLoc = [pathY(i+1),pathX(i+1)]; 
        end
        [closestIndex,closestX, closestY] = ComputeShortestPath(currentAdvLoc, nextAdvLoc, paths, flowsData);
        closestPaths{i} = paths{closestIndex};        
        EstimatedClosestAdvLoc(i,:)=[closestX, closestY];
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





