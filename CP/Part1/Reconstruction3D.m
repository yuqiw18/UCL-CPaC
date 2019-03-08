% Function for Synthetic Image 3D Reconstruction
%function result = Reconstruction3D(uvPatternSequence)

clearvars -except uvPatternSequence;

% Setup & Initialisation
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

%% 1. Light Patterns Decoding
% U: First 20 images 
% V: Remaining 20 images
decodedU = DecodeUV(uvPatternSequence(:,:,1:20));
decodedV = DecodeUV(uvPatternSequence(:,:,21:40));

%% 2. Unreliable Pixel Elimination
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
  

%% 4. Depth Map Computation

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

%% 5. Point Cloud Visualisation






%end


function decodedUV = DecodeUV(uvPatternSequence)
disp('@Decoding UV Pattern');
tic
    [height,width,~]=size(uvPatternSequence);    
    decodedUV = zeros(height,width,2);
    unknownThreshold = 0.05;
   
    % For each pixel
    for h=1:height
        for w=1:width
            
            % Decode using image difference
            codeWord = seq(h,w,1:2:40) - seq(h,w,2:2:40);
            
            codeWordU = codeWord(1:10);
            codeWordV = codeWord(11:20);
            
            % Convert to decimal values
            decodedUV(h,w,:) = [sum(codeWordU.*bin_seq'),sum(codeWordV.*bin_seq')];

            % Detect Unknown values
             if(sum(codeWord) == 0 && ~own_data)
                decodedUV(h,w,:) = -1; % Mark as unknown
             elseif(abs(sum(codeWord)) < T && own_data)
                decodedUV(h,w,:) = -1; % Mark as unknown
             end 
        end
    end
toc
end