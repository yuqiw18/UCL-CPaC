% Function for generating slow motion between two frames 
% Algorithm based on section 3 from this paper:
% Jiang, Huaizu & Sun, Deqing & Jampani, Varun & Yang, Ming-Hsuan & Learned-Miller, Erik & Kautz, Jan. (2017). 
% Super SloMo: High Quality Estimation of Multiple Intermediate Frames for Video Interpolation. 
% Retrieved from https://arxiv.org/pdf/1712.00080.pdf
% For detailed demonstration, please refer to the report.
function interpolatedFrames = MotionInterpolation(currentFrame, nextFrame, flow, frame)

    [width,height,channel] = size(currentFrame);
    interpolatedFrames = zeros(width,height,channel,frame);
    
    for f = 1:frame
        
        % Synthesize the intermediate image
        % The basic formula for interpolating two images is given by:
        % alpha: Contribution of the two selected images
        % f0: privious interpolated frame from warping with optical flow
        % f1: next interpolated frame from warping with optical flow
%         f0 = imwarp(currentFrame,-flow);  
%         f1 = imwarp(nextFrame, flow);  
%         interpolatedFrame = alpha*f0 + (1-alpha)*f1;
%         OR
%         interpolatedFrame = (1-alpha)*f0 + alpha*f1;
        
        % However, we need to interpolate n images between two given images,
        % Shift contribution to flow to reduce the artifact:        
        % f0 will have higher contribution in previous several frames but
        % will decrease and the contribution of f1 will increase so that
        % the artifact will be minimized.
        alpha = f/(frame+1);
        f0 = imwarp(currentFrame, alpha*(-flow));  
        f1 = imwarp(nextFrame, (1-alpha)*(flow));
        interpolatedFrame = (1-alpha)*f0 + alpha*f1;     
        
        % Save to output the image sequence
        interpolatedFrames(:,:,:,f) = interpolatedFrame;
    end
end