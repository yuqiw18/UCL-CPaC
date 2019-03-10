% Function for loading specific calibration matrix profile
function [cameraIntrinsic,cameraExtrinsic,cameraProjection,projectorIntrinsic,projectorExtrinsic,projectorProjection] = LoadCalibMatrixProfile(calibrationMatrix)
    switch calibrationMatrix
        case 'provided_synthetic'
            disp('*Using Provided Synthetic Calibration Matrix');
            cameraIntrinsic = [4786.25390625,0.00000000,1541.16491699;
                 0.00000000,4789.81884766,1036.94421387;
                 0.00000000,0.00000000,1.00000000];
             cameraExtrinsic = [0.99998617,-0.00475739,0.00223672,-0.10157115;
                 0.00382316,0.95016861,0.31171292,-0.10139455;
                 -0.00360820,-0.31170005,0.95017368,0.49625999];
             cameraProjection = [4780.62646484,-503.15124512,1475.07995605,278.67315674;
                  14.57075500,4227.92041016,2478.32568359,28.93240356;
                  -0.00360820,-0.31170005,0.95017368,0.49625999];
              projectorIntrinsic = [3680.39404297,0.00000000,591.75494385;
                 0.00000000,3672.32153320,393.62173462;
                 0.00000000,0.00000000,1.00000000];
             projectorExtrinsic = [0.72119248,0.44233182,-0.53312665,-0.14915472;
                 -0.36164442,0.89680630,0.25485638,-0.06104425;
                 0.59084243,0.00900178,0.80673677,1.36014771];
             projectorProjection = [3003.90649414,1633.28234863,-1484.72570801,255.92596436;
                  -1095.50610352,3296.90429688,1253.46362305,311.20962524;
                  0.59084243,0.00900178,0.80673677,1.36014771];              
        case 'calib_synthetic'
            disp('*Using Estimated Synthetic Calibration Matrix');
            cameraIntrinsic = [7575.74257718199,0,950.462527925728;0,5557.33651048385,-74.2905886684759;0,0,1];
            cameraExtrinsic = [-0.451192018497089,0.519616823513591,-0.725550907356726,905.617391173608;
                0.876976999418753,0.107533525115184,-0.468345901515944,1385.23317486611;
                -0.165339362932735,-0.847605390319072,-0.504210270985281,10587.4114975007];
            cameraProjection = cameraIntrinsic * cameraExtrinsic;
            
            projectorIntrinsic = [3760.16207644034,0,569.211395048940;0,3870.49298268604,534.297914776516;0,0,1];
            projectorExtrinsic = [0.371325341364319,0.928242164328734,0.0219994368787164,-1300.92414392579;
                0.873882817417158,-0.341377148756806,-0.346107589818144,-938.863043219109;
                -0.313761553227484,0.147743448819167,-0.937936864105093,12298.8747885383];
            projectorProjection = projectorIntrinsic*projectorExtrinsic;
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
end