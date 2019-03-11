%close all;
clearvars -except uvPatternSequence decodedUV;

% Setup & Initialisation
% Options: cube_T1 monkey_T1 notebook_T1 red_T1 sphere_T1 tablet_T1
% real_crayon_dalek real_tea 
% capture
[path, filename, prefix, first, last, digits, suffix, outputPath, threshold] = LoadImageSequenceProfile('capture');

% Load selected image sequence
if (exist('uvPatternSequence','var') == 0)
    disp("@Loading UV Patterns");
    tic  
        uvPatternSequence = load_sequence(strcat(path, filename, '/'), prefix, first, last, digits, suffix);
        uvPatternSequence = im2double(uvPatternSequence);
    toc
else
    disp("*UV Patterns Already Loaded");
end

%% 1,2. Light Patterns Decoding & Unreliable Pixel Elimination
if (exist('decodedUV','var') == 0)
    decodedUV = DecodeUV(uvPatternSequence, threshold);
else
    disp("*UV Already Decoded");
end

u = decodedUV(:,:,1);
v = decodedUV(:,:,2);
figure;
subplot(1,2,1),imagesc(u),title('U');
subplot(1,2,2),imagesc(v),title('V');

%% 3,4 Calibration Matrix Setup & Depth Map Computation
% Options: provided_synthetic calib_synthetic calib_real calib_capture
[depthMap,pointCloud] = ComputeDepthMap(decodedUV, 'calib_capture');
figure;
imagesc(depthMap),title('Depth Map');

% %% 5. Point Cloud Visualisation
SavePLY(pointCloud, strcat(outputPath,filename, '.ply'));
disp('>>Task Complete')