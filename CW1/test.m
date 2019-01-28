clc;
clear;

image = imread("footage/footage_625.png");
image2 = double(image);


[vertical,horizontal,sequenceLength]=size(image);


for v = 1:vertical   
    image2(v,:) = medfilt1(image2(v,:), 5);
end

imge2 = imsharpen(image2);

image2 = uint8(image2);




figure;
subplot(1,2,1);
imhist(image);


subplot(1,2,2);
imhist(image2);
