function output = labs3(path, prefix, first, last, digits, suffix)

% footage_001.png ~ footage_657.png

%
% Read a sequence of images and correct the film defects. This is the file 
% you have to fill for the coursework. Do not change the function 
% declaration, keep this skeleton. You are advised to create subfunctions.
% 
% Arguments:
%
% path: path of the files
% prefix: prefix of the filename
% first: first frame
% last: last frame
% digits: number of digits of the frame number
% suffix: suffix of the filename
%
% This should generate corrected images named [path]/corrected_[prefix][number].png
%
% Example:
%
% mov = labs3('../images','myimage', 0, 10, 4, 'png')
%   -> that will load and correct images from '../images/myimage0000.png' to '../images/myimage0010.png'
%   -> and export '../images/corrected_myimage0000.png' to '../images/corrected_myimage0010.png'
%

% Your code here
%% Image Initialisation

matrix = load_sequence(path, prefix, first, last, digits, suffix);


%% Scene Cut Detection






%% Global Flicker Correction






%% Blotch Correction






%% Vertical Artifacts Correction






%% Save the result
%save_sequence(matrix, path, prefix, first, digits);

output = 0;

end
