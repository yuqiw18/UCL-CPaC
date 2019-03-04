% Motion Interpolation 
% Algorithm based on this paper:
% Jiang, Huaizu & Sun, Deqing & Jampani, Varun & Yang, Ming-Hsuan & Learned-Miller, Erik & Kautz, Jan. (2017). 
% Super SloMo: High Quality Estimation of Multiple Intermediate Frames for Video Interpolation. 
% Retrieved from https://arxiv.org/pdf/1712.00080.pdf
% Formula (1)
function interpolatedFrames = MotionInterpolation(currentFrame, nextFrame, flow, frame)

    [width,height,channel] = size(currentFrame);
    interpolatedFrames = zeros(width,height,channel,frame);
    
    for f = 1:frame
        alpha = f/(frame+1);
        
        currentMotionCompensation = imwarp(currentFrame, alpha*(-flow));
        nextMotionCompensation = imwarp(nextFrame, (1-alpha)*(flow));
        
        interpolatedFrame = (1-alpha)*currentMotionCompensation + (alpha)*nextMotionCompensation;
        interpolatedFrames(:,:,:,f) = interpolatedFrame;
    end
end