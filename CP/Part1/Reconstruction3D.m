% Function for Synthetic Image 3D Reconstruction
%function result = Reconstruction3D(uvPatternSequence)

clearvars -except uvPatternSequence;

path = 'data/synthetic_data/';
prefix = 'cube_T1/';
first = 0;
last = 39;
digits = 4;
suffix = 'png';
outputPath = 'output';

if (exist('uvPatternSequence','var') == 0)
    disp("@Loading UV Patterns");
    tic  
        uvPatternSequence = load_sequence(path, prefix, first, last, digits, suffix);
        uvPatternSequence = double(uvPatternSequence);
    toc
else
    disp("@UV Patterns Already Loaded");
end

% U: First 20 images 
% V: Remaining 20 images
decodedU = DecodeUV(uvPatternSequence(:,:,1:20));
decodedV = DecodeUV(uvPatternSequence(:,:,21:40));

% generate mask to eliminate errors
mask_sum = decodedU(:,:,2)+decodedV(:,:,2);
mask = mask_sum > 5;
mask_3d = repmat(mask,[1,1,2]);

% u=uv_code(:,:,1) v=iv_code(:,:,2)
uv_code = cat(3,decodedU(:,:,1),decodedV(:,:,1));

% use mask to remove background
[height,width,~] = size(uvPatternSequence);
coded_pix = mask_3d .* uv_code + (1-mask_3d) .* zeros(height,width,2);

u = coded_pix(:,:,1);
v = coded_pix(:,:,2);

figure;
subplot(1,2,1),imagesc(u),title('u');
subplot(1,2,2),imagesc(v),title('v');
    
%     fprintf('Computeing 3D point cloud...');
%     [depth_map,point_cloud] = get_depth(coded_pix);
%     
%     mask = depth_map == -1;
%     depth_map = mask.*(max(depth_map(:)))+(1-mask).*depth_map;
%    
%     figure;
%     imagesc(depth_map),title('depth map');
%     %saveas(gcf, ['depth_map_',folder,'.jpg']);
%     fprintf('done\n');
%     
%     % save as PLY file -------------------------    
%     fprintf('Saving as PLY file...');
%     saveas_ply(point_cloud,folder);
%     fprintf('done\n');

%end


function decodedUV = DecodeUV(uvPatternSequence)
disp('@Decoding UV Pattern');
tic
    [height,width,~]=size(uvPatternSequence);    
    decodedUV = zeros(height,width,2);
   
    for h=1:height
        disp(h/height);
        for w=1:width
            disp('');
            % current_code represents binary 1-10 order
            current_code = zeros(1,10);
            
            diff_sum = 0;
            
            idx = 1;
            
            for k=1:2:20
                
                int1 = uvPatternSequence(h,w,k);     
                int2 = uvPatternSequence(h,w,k+1);
                
                % if white is before black,set as 1; otherwise, set as 0
                if int1<=int2
                    current_code(idx) = 0;
                else
                    current_code(idx) = 1;
                end
                
                idx = idx+1;
                diff_sum = diff_sum + abs(int1-int2);
            end
            
            % convert binary to decimal
            decodedUV(h,w,1) = bi2de(current_code)+1;
            decodedUV(h,w,2) = diff_sum;
        end
    end
toc
end