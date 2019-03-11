% Function for Decoding Light Patterns
function decodedUV = DecodeUV(uvPatternSequence, threshold)
disp('@Decoding UV Pattern');
tic  
    % Intialisation
    [height,width,~]=size(uvPatternSequence);    
    decodedUV = zeros(height,width,2);
    binary= [1 2 4 8 16 32 64 128 256 512];
    sumThreshold = threshold * 10;
    
    % Compute image difference between each pair of images
    imageDifference = uvPatternSequence(:,:,1:2:40) - uvPatternSequence(:,:,2:2:40);
    
    % First 10 results represent U, next 10 results represent V
    differenceU = imageDifference(:,:,1:10);
    differenceV = imageDifference(:,:,11:20);
    
    % Initialise storage for binary code words
    codewordU = zeros();
    codewordV = zeros();
    
    for h=1:height
        for w=1:width

            sumU = 0;
            sumV = 0;
            
            % Compute binary code word
            for d = 1:10

               if (differenceU(h,w,d)>threshold)
                   codewordU(d) = 1;
               elseif (differenceU(h,w,d)<-threshold)
                   codewordU(d) = 0;
               else
                   codewordU(d) = 0;
               end
               
               sumU = sumU + abs(differenceU(h,w,d));
               
               if (differenceV(h,w,d)>threshold)
                   codewordV(d) = 1;
               elseif (differenceV(h,w,d)<-threshold)
                   codewordV(d) = 0;
               else
                   codewordV(d) = 0;
               end
               
               sumV = sumV + abs(differenceV(h,w,d));
            end
            
            % Reject unreliable piexels & Convert to decimal values
            if (sumU < sumThreshold)
                decodedUV(h,w,1) = -1;
            else
                decodedUV(h,w,1) = sum(codewordU.*binary);
            end
                 
            if (sumV < sumThreshold)
                decodedUV(h,w,2) = -1;
            else    
                decodedUV(h,w,2) = sum(codewordV.*binary);
            end
            
        end
        
    end
toc 
end