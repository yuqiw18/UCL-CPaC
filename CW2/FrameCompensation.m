% function used for slow motion
% input ********************************************
% start_f -- start frame
% end_f -- end frame
% flow -- flow matrix for start and end frames
% f_num -- number of interpolated frames between two given frames
% Output *******************************************
% path_imgs -- interpolated sequence

function [path_imgs] = slow_motion(start_f, end_f, flow, f_num)
    [height, width, dim] = size(start_f);
    path_imgs = zeros(height, width, dim, f_num);
    step = 1/f_num;
    for f=1:f_num
        n=f*step;
        path_imgs(:,:,:,f) = get_interpolated_frame(start_f, end_f, flow, n);
    end
end
