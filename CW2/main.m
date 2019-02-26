clear;
clc;
warning('off','all');

path = 'data';
prefix = 'gjbLookAtTarget_';
first = 1;
last = 11;
digits = 4;
suffix = 'jpg';

outputPath = 'output';

% Load image sequence
rawImageSequence = load_sequence_color(path, prefix, first, last, digits, suffix);

% Downsize image size to 30% to match optical flow used in later stage
imageSequence = imresize(rawImageSequence, 0.3);

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

rawMatrix = distanceMatrix;

% 2. Convert the Distance Matrix into a graph.
for i = 1:w-1   
    average = sum(rawMatrix(i,i:w),'all')/(w-i); 
    
    disp(average);
    for j = i:w
        if (rawMatrix(i,j)>average)
            rawMatrix(i,j)=0;
        end
%         if (distanceMatrix(j,i)>average)
%             distanceMatrix(j,i)=0;
%         end
        
 rawMatrix(j,i)=0;

    end
end

sparseDistanceMatrix = sparse(rawMatrix);
graph = biograph(sparseDistanceMatrix, [],'ShowArrows','off','ShowWeights','on','LayoutType','equilibrium');
view(graph);

% 3. For each node in the graph, compute the shortest path.





% 4. Compute the advected location of point s in every other node, using
% optical flow.





% 5. Pick the path in Paths whose advected location comes closest to the selected point. 
% Render this path as the output image sequence for this user-drawn segment.




% 
% outputImageSequence = 
% implay(outputImageSequence)
