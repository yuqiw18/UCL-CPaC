%close all;
clearvars -except uvPatternSequence decodedUV;

% Setup & Initialisation
% path = 'data/real_data/';
% filename = 'real_crayon_dalek';
% prefix = 'IMG_';
% first = 9418;
% last = 9457;
% digits = 4;
% suffix = 'jpg';
% outputPath = 'output/';

%% Initialisation
% Image sequence parameters
path = 'data/synthetic_data/';
filename = 'cube_T1';
prefix = '';
first = 0;
last = 39;
digits = 4;
suffix = 'png';
outputPath = 'output/';

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

u = decodedUV(:,:,1);
v = decodedUV(:,:,2);
figure;
subplot(1,2,1),imagesc(u),title('U');
subplot(1,2,2),imagesc(v),title('V');

%% 3,4 Calibration Matrix Setup & Depth Map Computation
depthMap = ComputeDepthMap(decodedUV, 'provided_synthetic');

%% 5. Point Cloud Visualisation
%SavePLY(depthMap, strcat(outputPath,filename, '.ply'));
%pcwrite(depthMap,strcat(outputPath,filename),'PLYFormat','binary');

disp('>>Task Complete')
%% Core Functions

function depthMap = ComputeDepthMap(decodedUV, calibrationMatrix)
disp('@Computing Depth Map')
tic
switch calibrationMatrix
    case 'provided_synthetic'
        disp('*Using Provided Synthetic Calibration Matrix');
        SyntheticMatrixProfile;
    case 'calib_synthetic'
        disp('*Using Estimated Synthetic Calibration Matrix');
        CalibSyntheticMatrixProfile;
    case 'calib_real'
        disp('*Using Estimated Real Calibration Matrix');
        CalibRealMatrixProfile;
    case 'calib_capture'
        disp('*Using Estimated Captured Calibration Matrix');
        CalibCaptureMatrixProfile;
    otherwise
        disp('ERROR: Undefined Calibration Type');
        cameraIntrinsic = zeros(3,3);
        cameraExtrinsic = zeros(3,4);
        cameraProjection = zeros(3,4);     
        projectorIntrinsic = zeros(3,3);
        projectorExtrinsic = zeros(3,4);
        projectorProjection = zeros(3,4);     
end             
%     cameraIntrinsic = [[4786.25390625,0.00000000,1541.16491699]
%                  [0.00000000,4789.81884766,1036.94421387]
%                  [0.00000000,0.00000000,1.00000000]];
%              
%     cameraExtrinsic =[[0.954846897455116,-0.243793610738522,0.339600222517088,0.103748],
%                  [0.295599497778259,0.722241740754021,-1.25058035299751,0.250543],
%                  [0.0298050219523619,0.647249212649715,1.52339110853483,-0.439700]];
% 
%     cameraProjection = cameraIntrinsic*cameraExtrinsic;
%     
%     projectorIntrinsic =[[3680.39404297,0.00000000,591.75494385],
%                  [0.00000000,3672.32153320,393.62173462],
%                  [0.00000000,0.00000000,1.00000000]];
%              
%     projectorExtrinsic =    [[-0.389657586760984,-0.163708660668265,1.81258537950830,-0.718140],
%                  [-0.852751906358342,-0.307524981505965,-0.844375679308612,0.108477],
%                  [0.347823448435870,-0.937351513666009,-0.0395468606198075,-1.161242]];
%              
%     projectorProjection = projectorIntrinsic*projectorExtrinsic;
              
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