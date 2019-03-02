%%
%
path = 'mydata';
prefix = 'mydata_';
first = 0;
last = 24;
digits = 4;
suffix = 'png';
outputPath = 'output';

%
flowImageSequence = load_sequence_color(path, prefix, first, last, digits, suffix);

[height,width,~,N] = size(flowImageSequence);
flows_a = zeros(height, width, 2, N*(N-1)/2);

% set parameters
alpha = 0.01;     % the regularization weight
ratio = 0.75;      % the downsample ratio
minWidth = 32;     % the width of the coarest level
nOuterFPIterations = 8;     % the number of outer fixed point iterations
nInnerFPIterations = 1;     % the number of inner fixed point iterations
nSORIterations = 32;        % the number of SOR iteration

parameters = [alpha, ratio, minWidth, nOuterFPIterations, nInnerFPIterations, nSORIterations];

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

save('flow/myflows.mat','flows_a');