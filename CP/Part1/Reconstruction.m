%close all;
clearvars -except uvPatternSequence decodedUV;

% Setup & Initialisation
[path, filename, prefix, first, last, digits, suffix, outputPath] = LoadImageSequenceProfile('cube_T1');

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
    decodedUV = DecodeUV(uvPatternSequence);
else
    disp("*UV Already Decoded");
end

% u = decodedUV(:,:,1);
% v = decodedUV(:,:,2);
% figure;
% subplot(1,2,1),imagesc(u),title('U');
% subplot(1,2,2),imagesc(v),title('V');

%% 3,4 Calibration Matrix Setup & Depth Map Computation
depthMap = ComputeDepthMap(decodedUV, 'calib_synthetic');

%% 5. Point Cloud Visualisation
SavePLY(depthMap, strcat(outputPath,filename, '.ply'));

disp('>>Task Complete')