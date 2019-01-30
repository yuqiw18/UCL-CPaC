function output = labs3(path, prefix, first, last, digits, suffix)

% footage_001.png ~ footage_657.png

%
% Read a sequence of images and correct the film defects. This is the file 
% you have to fill for the coursework. Do not change the function 
% declaration, keep this skeleton. You are advised to create subfunctions.
% 
% Arguments:
%
% path: path of the files
% prefix: prefix of the filename
% first: first frame
% last: last frame
% digits: number of digits of the frame number
% suffix: suffix of the filename
%
% This should generate corrected images named [path]/corrected_[prefix][number].png
%
% Example:
%
% mov = labs3('../images','myimage', 0, 10, 4, 'png')
%   -> that will load and correct images from '../images/myimage0000.png' to '../images/myimage0010.png'
%   -> and export '../images/corrected_myimage0000.png' to '../images/corrected_myimage0010.png'
%

% Your code here
%% Image Initialisation

rawImageSequence = load_sequence(path, prefix, first, last, digits, suffix);
[vertical,horizontal,sequenceLength]=size(rawImageSequence);
imageSequence = im2double(rawImageSequence);

imageSequenceEdge = zeros(size(imageSequence));
motionMaskSequence = zeros(size(imageSequence));
blotchMaskSequence = zeros(size(imageSequence));

sharpenFilter = [0,0,0; 0,2,0; 0,0,0];
meanFilter = [1/9,1/9,1/9;1/9,1/9,1/9;1/9,1/9,1/9];
laplacianFilter = fspecial('laplacian',0);

sceneCutThreshold = 45000;

%% Scene Cut Detection

[sceneCutFrames, sceneClipPoints] = SceneCutDetection(imageSequence,first,last,sceneCutThreshold);

imageSequence = OverlayText(imageSequence,sceneCutFrames);

%% Global Flicker Correction






%% Blotch Correction






%% Vertical Artifacts Correction






%% Save the result
%save_sequence(matrix, path, prefix, first, digits);

output = 0;

end

function [sceneCutFrames, sceneClipPoints] = SceneCutDetection(imageSequence,first,last,sceneCutThreshold)

sceneCutFrames = [];

sequenceLength = length(imageSequence);

for i = 1+1 : sequenceLength-1
    %if (i ~= 1 && i+1<sequenceLength)       
                      
        %
        previousFrameDiff = (sum(abs(imageSequence(:,:,i)-imageSequence(:,:,i-1)),'all'));
        nextFrameDiff = (sum(abs(imageSequence(:,:,i)-imageSequence(:,:,i+1)),'all'));          
        transitionDiff = abs(previousFrameDiff - nextFrameDiff);
              
        if (transitionDiff > sceneCutThreshold)     
            sceneCutFrames = [sceneCutFrames, i+first-1];
            %disp(strcat('Scene Cut Detected @', int2str(i+first-1)));
        end
end

sceneClipPoints = sceneCutFrames;
% Append first frame
sceneClipPoints = [first, sceneClipPoints];
% Append last frame
sceneClipPoints = [sceneClipPoints, last];
% Remove "redundant" cut scene frame
sceneCutFrames = sceneCutFrames(2:2:length(sceneCutFrames));

disp('Scene Cut Detected @');
disp(sceneCutFrames);

%sceneClipPointsCount = length(sceneClipPoints);

end

function imageSequence = OverlayText(imageSequence,sceneCutFrames)
    for f = 1 : length(sceneCutFrames)
        detectedFrame = imageSequence(:,:,sceneCutFrames(f));
        detectedFrame = insertText(detectedFrame,[0,0],'Scene Cut Detected','FontSize',24); 
        imageSequence(:,:,sceneCutFrames(f)) = rgb2gray(detectedFrame);
    end
end


