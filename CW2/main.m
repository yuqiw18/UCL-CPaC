clear;
clc;
warning('off','all');

% Setup
path = 'data';
prefix = 'gjbLookAtTarget_';
first = 1;
last = 71;
digits = 4;
suffix = 'jpg';
outputPath = 'output';

% Initialisation
flowsFile = load('flow/flows.mat');
flowsData = flowsFile.flows_a;
selectedImageIndex = 5;

% Load image sequence
rawImageSequence = load_sequence_color(path, prefix, first, last, digits, suffix);

% Downsize image size to 30% to match optical flow used in later stage
imageSequence = imresize(rawImageSequence, 0.3);
imageCount = size(imageSequence,4);

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

sparseDistanceMatrix = sparse(connectionMatrix);
% graph = biograph(sparseDistanceMatrix, [],'ShowArrows','off','ShowWeights','on','LayoutType','equilibrium');

% graph = graphminspantree(sparseDistanceMatrix);
% view(graph);

% 3. For each node in the graph, compute the shortest path.
[~,paths,~] = graphshortestpath(sparseDistanceMatrix,selectedImageIndex);
% disp(paths);

% 4. Compute the advected location of point s in every other node, using
% optical flow.

% [height, width, pair, ~] = size(flowsData);
% flowValueList = zeros(imageCount,height,width,pair);
% 
% for i=1:size(paths,2)
%     
%    currentPath = paths(i);
%    
%    currentPath = currentPath{1};
%    
%    totalFlow = zeros(height,width,pair);
%    
%    for f = 1:size(currentPath,2)
%        currentFrameIndex = currentPath(f);
%        
%        if (f+1 > size(currentPath,2))
%            nextFrameIndex = currentPath(f);
%        else
%            nextFrameIndex = currentPath(f+1);
%        end
%        
%        if (currentFrameIndex > nextFrameIndex )
%            
%            k = (currentFrameIndex-1)*(currentFrameIndex-2)/2 + nextFrameIndex;        
%            flow = flowsData(:,:,:,k);
%            totalFlow = totalFlow + flow;
%        elseif (currentFrameIndex < nextFrameIndex)
%            k = (nextFrameIndex-1)*(nextFrameIndex-2)/2+currentFrameIndex;
%            flow = -flowsData(:,:,:,k);
%            totalFlow = totalFlow + flow;
%        else
%            flow = zeros(height,width,pair);
%            totalFlow = totalFlow + flow;
%        end 
%    end
%    flowValueList(i,:,:,:) = totalFlow;
% end

% 5. Pick the path in Paths whose advected location comes closest to the selected point. 
% Render this path as the output image sequence for this user-drawn segment.

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

outputIndex = ConvertPathToIndex(closestPaths);

outputImageSequence = zeros(360, 480, 3, length(outputIndex));

for f = 1:length(outputIndex)
   outputImageSequence(:,:,:,f) = imageSequence(:,:,:,outputIndex(f)); 
end

implay(outputImageSequence);

% outputImageSequence = 
% implay(outputImageSequence)
