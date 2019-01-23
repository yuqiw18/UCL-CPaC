clear;
clc;

%
% cd /Users/yuqi/Desktop/UCL-CPaC/CW1/footage_raw 
% ffmpeg -i footage_%03d.png video.mp4

%
path = 'footage';
prefix = 'footage_';
first = 1;
last = 657;
digits = 3;
suffix = 'png';

% footage_001.png ~ footage_657.png

% Your code here
%% Image Initialisation

rawImageSequence = load_sequence(path, prefix, first, last, digits, suffix);

[~,~,sequenceLength]=size(rawImageSequence);

imageSequence = im2double(rawImageSequence);

%% Scene Cut Detection
threshold = 45000;
for i = 1 : sequenceLength 
    if (i ~= 1 && i+1<sequenceLength)       
                      
        %
        previousFrameDiff = (sum(abs(imageSequence(:,:,i)-imageSequence(:,:,i-1)),'all'));
        nextFrameDiff = (sum(abs(imageSequence(:,:,i)-imageSequence(:,:,i+1)),'all'));
            
        transitionDiff = abs(previousFrameDiff - nextFrameDiff);
              
        if (transitionDiff > threshold)     
            disp(strcat('Scene Cut Detected @', int2str(i+first-1)));
        end
        
        % Debug Info
%         disp(strcat("Previous Frame Diff: ",num2str(abs(previousFrameDiff))));
%         disp(strcat("Next Frame Diff: ",num2str(abs(nextFrameDiff))));
%         disp(strcat("Transition Diff: ", num2str(transitionDiff)));
%         disp(" ");        
    end 
end

%% Global Flicker Correction






%% Blotch Correction






%% Vertical Artifacts Correction







%% Camera Shake Calibration







%% Save the result
%save_sequence(matrix, path, prefix, first, digits);