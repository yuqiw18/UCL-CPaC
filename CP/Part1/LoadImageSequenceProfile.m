% Function for loading image sequence parameters
function [path, filename, prefix, first, last, digits, suffix, outputPath, threshold] = LoadImageSequenceProfile(datasetName)
switch datasetName
    case 'cube_T1'
    path = 'data/synthetic_data/';
    filename = 'cube_T1';
    prefix = '';
    first = 0;
    last = 39;
    digits = 4;
    suffix = 'png';
    outputPath = 'output/';
    threshold = 0.2;
    
    case 'monkey_T1'
    path = 'data/synthetic_data/';
    filename = 'monkey_T1';
    prefix = '';
    first = 0;
    last = 39;
    digits = 4;
    suffix = 'png';
    outputPath = 'output/';  
    threshold = 0.16;
        
    case 'notebook_T1'
    path = 'data/synthetic_data/';
    filename = 'notebook_T1';
    prefix = '';
    first = 0;
    last = 39;
    digits = 4;
    suffix = 'png';
    outputPath = 'output/'; 
    threshold = 0.16;
        
    case 'red_T1'
    path = 'data/synthetic_data/';
    filename = 'red_T1';
    prefix = '';
    first = 0;
    last = 39;
    digits = 4;
    suffix = 'png';
    outputPath = 'output/';
    threshold = 0.16;
            
    case 'sphere_T1'
    path = 'data/synthetic_data/';
    filename = 'sphere_T1';
    prefix = '';
    first = 0;
    last = 39;
    digits = 4;
    suffix = 'png';
    outputPath = 'output/';    
    threshold = 0.16;
    
    case 'tablet_T1'
    path = 'data/synthetic_data/';
    filename = 'tablet_T1';
    prefix = '';
    first = 0;
    last = 39;
    digits = 4;
    suffix = 'png';
    outputPath = 'output/';
    threshold = 0.16;
    
    case 'real_crayon_dalek'
    path = 'data/real_data/';
    filename = 'real_crayon_dalek';
    prefix = 'IMG_';
    first = 9418;
    last = 9457;
    digits = 4;
    suffix = 'jpg';
    outputPath = 'output/';
    threshold = 0.16;
    
    case 'real_tea'
    path = 'data/real_data/';
    filename = 'real_tea';
    prefix = 'IMG_';
    first = 9377;
    last = 9416;
    digits = 4;
    suffix = 'jpg';
    outputPath = 'output/';
    threshold = 0.16;
        
    case 'capture'
    path = 'data/capture_data/';
    filename = 'head_coffee';
    prefix = 'IMG_';
    first = 3296;
    last = 3335;
    digits = 4;
    suffix = 'jpg';
    outputPath = 'output/';
    threshold = 0.08;
           
    otherwise
    path = 'data/synthetic_data/';
    filename = 'cube_T1';
    prefix = '';
    first = 0;
    last = 39;
    digits = 4;
    suffix = 'png';
    outputPath = 'output/';
    threshold = 0.16;
end
end