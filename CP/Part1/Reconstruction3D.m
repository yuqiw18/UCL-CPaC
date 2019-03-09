% Function for Synthetic Image 3D Reconstruction

close all;
clearvars -except uvPatternSequence decodedUV;

% Setup & Initialisation
path = 'data/synthetic_data/';
filename = 'cube_T1';
prefix = '';
first = 0;
last = 39;
digits = 4;
suffix = 'png';
outputPath = 'output';

if (exist('uvPatternSequence','var') == 0)
    disp("@Loading UV Patterns");
    tic  
        uvPatternSequence = load_sequence(strcat(path, filename, '/'), prefix, first, last, digits, suffix);
        uvPatternSequence = double(uvPatternSequence);
    toc
else
    disp("@UV Patterns Already Loaded");
end

%% 1,2. Light Patterns Decoding & Unreliable Pixel Elimination
if (exist('decodedUV','var') == 0)
    decodedUV = DecodeUV(uvPatternSequence);
else
    disp("@UV Already Decoded");
end

% u = decodedUV(:,:,1);
% v = decodedUV(:,:,2);
% figure;
% subplot(1,2,1),imagesc(u),title('U');
% subplot(1,2,2),imagesc(v),title('V');

%% 3,4 Calibration Matrix Setup & Depth Map Computation
depthMap = ComputeDepthMap(decodedUV, 0);

%% 5. Point Cloud Visualisation
%SavePLY(depthMap, strcat(filename, '.ply'));

%% Core Functions
function decodedUV = DecodeUV(uvPatternSequence)
disp('@Decoding UV Pattern');
tic
    [height,width,~]=size(uvPatternSequence);    
    decodedUV = zeros(height,width,2);
    binary= [1 2 4 8 16 32 64 128 256 512];
     
    codeWord = uvPatternSequence(:,:,1:2:40) - uvPatternSequence(:,:,2:2:40);
    codeWordU = codeWord(:,:,1:10);
    codeWordV = codeWord(:,:,11:20);
   
    for h=1:height
        %disp(h/height*100);
        for w=1:width         
            decodedUV(h,w,1) = sum(reshape(codeWordU(h,w,:)>=0,1,10).*binary);
            decodedUV(h,w,2) = sum(reshape(codeWordV(h,w,:)>=0,1,10).*binary);
            if sum(codeWord(h,w,:)) == 0
                decodedUV(h,w,:) = -1;
            end
        end
    end
toc 
end

function depthMap = ComputeDepthMap(decodedUV, calibrationMatrix)
disp('@Computing Depth Map')
tic
    cameraIntrinsic = [[   4786.25390625,       0.00000000,    1541.16491699]
                 [      0.00000000,    4789.81884766,    1036.94421387]
                 [      0.00000000,       0.00000000,       1.00000000]];
             
    cameraExtrinsic =[[      0.99998617,      -0.00475739,       0.00223672,      -0.10157115],
                 [      0.00382316,       0.95016861,       0.31171292,      -0.10139455],
                 [     -0.00360820,      -0.31170005,       0.95017368,       0.49625999]];
             
    cameraProjection = [[   4780.62646484,    -503.15124512,    1475.07995605,     278.67315674],
                  [     14.57075500,    4227.92041016,    2478.32568359,      28.93240356],
                  [     -0.00360820,      -0.31170005,       0.95017368,       0.49625999]];

    projectorIntrinsic =[[   3680.39404297,       0.00000000,     591.75494385],
                 [      0.00000000,    3672.32153320,     393.62173462],
                 [      0.00000000,       0.00000000,       1.00000000]];
             
    projectorExtrinsic =    [[      0.72119248,       0.44233182,      -0.53312665,      -0.14915472],
                 [     -0.36164442,       0.89680630,       0.25485638,      -0.06104425],
                 [      0.59084243,       0.00900178,       0.80673677,       1.36014771]];
             
    projectorProjection = [[   3003.90649414,    1633.28234863,   -1484.72570801,     255.92596436],
                  [  -1095.50610352,    3296.90429688,    1253.46362305,     311.20962524],
                  [      0.59084243,       0.00900178,       0.80673677,       1.36014771]];
         
    [w,h,~] = size(decodedUV);
    depthMap = zeros(w,h,3);

    % Load Intrinsic Matrix
    K1 = cameraIntrinsic;
    K2 = projectorIntrinsic;

    % Compute Rotation/Translation matrices
    R1 = cameraExtrinsic(1:3,1:3);
    T1 = cameraExtrinsic(1:3,4);
    R2 = projectorExtrinsic(1:3,1:3);
    T2 = projectorExtrinsic(1:3,4);

    % For each point solve the linear system
    for i=1:w
        %disp(i/w*100);
        for j=1:h
            if(decodedUV(i,j,1)~=-1)
                % Compute A and b
                A = zeros(4,3);
                b = zeros(4,1);

                % Convert points to normalized camera coordinates
                p1 = K1\[j;i;1];
                p2 = K2\[decodedUV(i,j,2); decodedUV(i,j,1);1];

                % Compute linear constraints -> Camera
                A(1,:) = [R1(3,1)*p1(1)-R1(1,1), R1(3,2)*p1(1)-R1(1,2), R1(3,3)*p1(1)-R1(1,3)];
                A(2,:) = [R1(3,1)*p1(2)-R1(2,1), R1(3,2)*p1(2)-R1(2,2), R1(3,3)*p1(2)-R1(2,3)];
                b(1:2) = [T1(1)-T1(3)*p1(1), T1(2)-T1(3)*p1(2)];

                % Compute linear constraints -> Projector (~ camera no.2)
                A(3,:) = [R2(3,1)*p2(1)-R2(1,1), R2(3,2)*p2(1)-R2(1,2), R2(3,3)*p2(1)-R2(1,3)];
                A(4,:) = [R2(3,1)*p2(2)-R2(2,1), R2(3,2)*p2(2)-R2(2,2), R2(3,3)*p2(2)-R2(2,3)];
                b(3:4) = [T2(1)-T2(3)*p2(1), T2(2)-T2(3)*p2(2)];

                % Compute the least squares solution
                x = A\b;

                % Store Depth
                depthMap(i,j,:) = cameraExtrinsic(1:3,1:4)*[x;1];
            end
        end
    end

toc
end

