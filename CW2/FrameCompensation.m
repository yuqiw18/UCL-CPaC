% function used to get interpolated frame for slow motion 
% Input *****************************************************
% start_f -- start frame
% end_f -- end frame
% flow -- corresponded flow matrix for start and end frames
% n -- nth step
% Output ****************************************************
% interpolated_f -- interpolated frame
function [interpolated_f] = get_interpolated_frame(start_f,end_f,flow,n)
   start_f = double(start_f);
   end_f = double(end_f);
   flow = double(flow);
   
   [height,width,dim] = size(start_f);
   [cord_x, cord_y] = meshgrid(1:width, 1:height);
   
   flow_ux = flow(:,:,1);
   flow_uy = flow(:,:,2);
   
   % get coordinate of original pixels after 'flow' from start frame
   cord_x_forward = cord_x + n*flow_ux;
   cord_y_forward = cord_y + n*flow_uy;
   
   % interpolate original coordinates from start frame
   flow_ux = griddata(cord_x_forward, cord_y_forward, flow_ux, cord_x, cord_y);
   flow_uy = griddata(cord_x_forward, cord_y_forward, flow_uy, cord_x, cord_y);
   
   % get coordinate of original pixels before 'flow' from end frame
   cord_x_back = cord_x + (1-n)*flow_ux;
   cord_y_back = cord_y + (1-n)*flow_uy;
   
   % interpolate from start frame (forward)
   start_f_forward = zeros(height, width, dim);
   
   % interpolate from end frame (back)
   end_f_back = zeros(height, width, dim);
   
   for d=1:dim
       start_f_forward(:,:,d) = griddata(cord_x_forward, cord_y_forward, start_f(:,:,d), cord_x, cord_y);
       end_f_back(:,:,d) = interp2(cord_x, cord_y, end_f(:,:,d), cord_x_back, cord_y_back, 'linear', NaN);
   end
   
   % replace NAN interpolated pixels
   nan_idx_s = find(isnan(start_f_forward));
   nan_idx_e = find(isnan(end_f_back));
   start_f_forward(nan_idx_s) = start_f(nan_idx_s);   
   end_f_back(nan_idx_e) = end_f(nan_idx_e);
   
   % blend frmaes from both two directions
   interpolated_f = (1-n)*start_f_forward + n*end_f_back;
end
