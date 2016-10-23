function mask = SubtractDominantMotion(image1, image2)
% 
% % input - image1 and image2 form the input image pair
% % output - mask is a binary image of the same size
% 
% M = LucasKanadeAffine(image1,image2);
% thresh = 70;
% tform = affine2d(M');
% warp_im1 = imwarp(double(image1),tform);
% result_im1 = warp_im1(1:size(image1,1),1:size(image1,2));
% error = abs(double(image2)-result_im1);
% 
% error(error<=thresh)=1;
% error(error>thresh)=0;
% small_removed = bwareaopen(error,50);
% SE = strel('square',5);
% erodedI = imerode(small_removed,SE);
% small_removed = bwareaopen(erodedI,100);
% 
% mask = small_removed;
% 
% 
% imshow(mask)

       It=image1;
        It1=image2;
        maskthreshold=30;

        M = LucasKanadeAffine(It, It1)
        t = affine2d(M');
        
        warpped_im = imwarp( It, t, 'OutputView' ,imref2d(size(It1)) );
      
        
        diff_image = warpped_im - It1;
        
        [ It_row , It_col ]=size( It ); 
        
        mask=zeros(It_row,It_col);        
        
        gretIndex = abs(diff_image) <  maskthreshold;
        lessIndex = abs(diff_image) >= maskthreshold;
        
        mask(gretIndex) = 0;        
        mask(lessIndex) = 1;
        se_mask = strel('disk',8);  
        mask = bwareaopen(mask,7);
        mask = imdilate(mask,se_mask);
        mask = imerode(mask,se_mask);
        
        