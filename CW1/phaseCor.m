clc;
clear;

referenceFrame = rgb2gray(im2single(imread('1.png')));
currentFrame = rgb2gray(im2single(imread('2.png')));

[row, col, ~] = size(referenceFrame);

windowRow = 2 * row -1;
windowCol = 2 * col -1;

Ga = fft2(referenceFrame,windowRow,windowCol);
Gb = fft2(currentFrame,windowRow,windowCol);

GbConj = conj(Gb);
R = (Ga .* GbConj)./abs(Ga .* GbConj); 
r = ifft2(R);

r = fftshift(r);

subpixel = true;
[xpeak,ypeak,peakVal] = findpeak(r, true);

peakValue = max(abs(r(:)));
[peakRowM,peakColM] = find(r==peakValue);

[peakRowW, peakColW] = findWeightedPeak(peakRowM, peakColM, r, 2);

if all(r(:) == peakValue)  
    frameShift = affine2d([1, 0, 0; 0, 1, 0; 0, 0, 1]);
else
 
    shiftRow = (1 + (windowRow-1)/2); 
    shiftCol = (1 + (windowCol-1)/2);

    frameShiftRow = peakRowW-shiftRow;
    frameShiftCol = peakColW-shiftCol;

    frameShift = affine2d([1, 0, 0; 0, 1, 0; frameShiftCol, frameShiftRow, 1]);
end

d = imregcorr(currentFrame, referenceFrame, 'translation'); 


function [weightedPeakRow, weightedPeakCol] = findWeightedPeak(peakRow, peakCol, r, span)

[rowLimit, colLimit] = size(r);

weight = 0;
weightedX = 0;
weightedY = 0;

startRow = peakRow - span;
endRow = peakRow + span;
startCol = peakCol - span;
endCol = peakCol + span;

if (peakRow - span < 0)
    startRow = 0;
end

if (peakRow + span > rowLimit)
    endRow = rowLimit;
end

if (peakCol - span < 0)
    startCol = 0;
end

if (peakCol + span > colLimit)
    endCol = colLimit;
end

for i=startRow:endRow
   for j = startCol:endCol
        weightedX = weightedX + j * r(i,j);
        weightedY = weightedY + i * r(i,j);
        weight = weight + r(i,j);
   end 
end

weightedPeakCol = weightedX/weight;
weightedPeakRow = weightedY/weight;

end
