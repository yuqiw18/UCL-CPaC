%function CameraCalibration(mode)

    calibrationMode = 'synthetic';

    switch calibrationMode
        case 'synthetic'
            for sn=0:2:4
            proj_name = ['data/synthetic_calibration/000',num2str(sn-1),'.png'];
            cam_name = ['data/synthetic_calibration/000',num2str(sn),'.png'];
            proj_img = imread(proj_name);
            cam_img = imread(cam_name);
            save_path = ['data/synthetic_calibration/reproject_000',num2str(sn-1),'.png'];

            [proj_img_out, cam_img_out] = reprojection(proj_img, cam_img, 1024, 768, 378, 277, 270, true);
            imwrite(cam_img_out, save_path);
            end
        case 'real'
            for rn=1:9
            proj_name = ['data/real_calibration/IMG_932',num2str(rn),'.jpg'];
            cam_name = ['data/real_calibration/IMG_932',num2str(rn),'.jpg'];
            proj_img = imread(proj_name);
            cam_img = imread(cam_name);
            save_path = ['data/real_calibration/reproject_932',num2str(rn),'.jpg'];
            
            [proj_img_out, cam_img_out] = reprojection(proj_img, cam_img, 1024, 768, 518, 120, 299, true);
            imwrite(cam_img_out,save_path);
            end
        case 'own'
            for on=1:6
            proj_name = ['data/own_calibration/IMG_28',num2str(on+87),'.jpg'];
            cam_name = ['data/own_calibration/IMG_28',num2str(on+87),'.jpg'];
            proj_img = imread(proj_name);
            cam_img = imread(cam_name);
            save_path = ['data/own_calibration/reproject_28',num2str(on+87),'.jpg'];
            
            [proj_img_out, cam_img_out] = reprojection(proj_img, cam_img, 1920, 1080, 705, 285, 512, true);
            imwrite(cam_img_out,save_path);
            end
        otherwise
            disp('ERROR: Undefined Calibration Mode');
    end
   
%% 1,2. Projection Matrix Estimation & Reprojection







%% 3. Point Cloud Generation & Visualisation
    






%end