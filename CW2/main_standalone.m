% clear;
% clc;
% warning('off','all');

% Setup
path = 'data';
prefix = 'gjbLookAtTarget_';
first = 1;
last = 71;
digits = 4;
suffix = 'jpg';
outputPath = 'output';

% Initialisation

disp("@Load Optical Flow Data");
tic  
    flowsFile = load('flow/flows.mat');
    flowsData = flowsFile.flows_a;
    selectedImageIndex = 5;
toc

%Load image sequence

disp("@Load Image Sequence");
tic  
    rawImageSequence = load_sequence_color(path, prefix, first, last, digits, suffix);
    imageSequence = imresize(rawImageSequence, 0.3);
toc

% Downsize image size to 30% to match optical flow used in later stage

imageCount = last - first + 1;
% imageCount = size(imageSequence,4);

%% Basic Section

% 1. Compute a Distance Matrix that encodes the similarity in appearance
% between all the frames in the collection.
distanceMatrix = ComputeDistanceMatrix(imageSequence);
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

imshow(imageSequence(:,:,:,selectedImageIndex)),title('Draw a path with at least 5 points');
pointCount = 0;
while(pointCount<5)
    
    [pathX, pathY]=getline();
    hold on;
    plot(pathX,pathY,'g*-');
    
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

outputIndex = ConvertPathToIndex(closestPaths);
outputImageSequence = zeros(360, 480, 3, length(outputIndex));

for f = 1:length(outputIndex)
   outputImageSequence(:,:,:,f) = imageSequence(:,:,:,outputIndex(f)); 
end

implay(outputImageSequence);





