% Generate optical flow from given images
%% Image parameters
path = 'mydata';
prefix = 'mydata_';
first = 0;
last = 24;
digits = 4;
suffix = 'png';
outputPath = 'output';

%% Initialisation
flowImageSequence = load_sequence_color(path, prefix, first, last, digits, suffix);
[height,width,~,N] = size(flowImageSequence);
flows_a = zeros(height, width, 2, N*(N-1)/2);

%% Optical flow parameters
% Regularization weight
alpha = 0.01;
% Downsample ratio
ratio = 0.75;
% Coarest level width
minWidth = 32;
% Outer fixed point iteration
nOuterFPIterations = 8;
% Inner fixed point iteration
nInnerFPIterations = 1;
% SOR iteration
nSORIterations = 32;   
parameters = [alpha, ratio, minWidth, nOuterFPIterations, nInnerFPIterations, nSORIterations];

%% Generate optical flow
k = 1;        
for i=2:N
   for j=1:i-1
       firstImage = flowImageSequence(:,:,:,i);
       secondImage = flowImageSequence(:,:,:,j);
       [x, y, ~] = Coarse2FineTwoFrames(firstImage, secondImage, parameters);
       flows_a(:,:,1,k) = x;
       flows_a(:,:,2,k) = y;
       k = k + 1;
    end
end 

%% Save data
save('flow/myflows.mat','flows_a');