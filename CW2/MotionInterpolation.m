function interpolatedFrames = MotionInterpolation(currentFrame, nextFrame, flow, frame)

    [width,height,channel] = size(currentFrame);
    interpolatedFrames = zeros(width,height,channel,frame);
    
    for f = 1:frame
        step = f/(frame+1);
        currentMotionCompensation = imwarp(currentFrame, step*(-flow));
        nextMotionCompensation = imwarp(nextFrame, (1-step)*(flow));
        interpolatedFrame = (1-step)*currentMotionCompensation + (step)*nextMotionCompensation;
        interpolatedFrames(:,:,:,f) = interpolatedFrame;
    end
end