% Function for Decoding Light Patterns
function decodedUV = DecodeUV(uvPatternSequence)
disp('@Decoding UV Pattern');
tic  
    % Intialisation
    [height,width,~]=size(uvPatternSequence);    
    decodedUV = zeros(height,width,2);
    binary= [1 2 4 8 16 32 64 128 256 512];
    threshold = 0.02;
    
    % Compute image difference between each pair
    codeWord = uvPatternSequence(:,:,1:2:40) - uvPatternSequence(:,:,2:2:40);
    
    % First 10 results represent U, next 10 results represent V
    codeWordU = codeWord(:,:,1:10);
    codeWordV = codeWord(:,:,11:20);
    
    % Initialise storage for binary code words
    currentU = zeros();
    currentV = zeros();
    
    for h=1:height
        for w=1:width
            
            isReliableU = true;
            isReliableV = true;
            sumU = 0;
            sumV = 0;
            
            % Compute binary code word
            for d = 1:10
               
               
               if (codeWordU(h,w,d)>threshold)
                   currentU(d) = 1;
               elseif (codeWordU(h,w,d)<-threshold)
                   currentU(d) = 0;
               else
                   currentU(d) = 0;
                   isReliableU = false;
               end  
               sumU = sumU+codeWordU(h,w,d);
               
               if (codeWordV(h,w,d)>threshold)
                   currentV(d) = 1;
               elseif (codeWordV(h,w,d)<-threshold)
                   currentV(d) = 0;
               else
                   currentV(d) = 0;
                   isReliableV = false;
               end
               sumV = sumV + codeWordV(h,w,d);
            end
            
            if (sumU > -threshold && sumU < threshold)
                decodedUV(h,w,1) = -1;
                disp('unrU')
            else
                decodedUV(h,w,1) = sum(currentU.*binary);
            end
                 
            if (sumV > -threshold && sumV < threshold)
                decodedUV(h,w,2) = -1;
                disp('unrV')
            else    
                decodedUV(h,w,2) = sum(currentV.*binary);
            end
            
            
%             % Reject unreliable piexels & Convert to decimal values
%             if (isReliableU)
%                 decodedUV(h,w,1) = sum(currentU.*binary);
%             else
%                 decodedUV(h,w,1) = -1;
%             end
%      
%             if (isReliableV)
%                 decodedUV(h,w,2) = sum(currentV.*binary);
%             else
%                 decodedUV(h,w,2) = -1;
%             end
        end
        
    end
toc 
end