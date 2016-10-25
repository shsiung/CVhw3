function mask = SubtractDominantMotion(image1, image2)
% 
% % input - image1 and image2 form the input image pair
% % output - mask is a binary image of the same size

thresh=30;

M = LucasKanadeAffine(image1, image2);
t = affine2d(M');

warpped_im = imwarp( image1, t, 'OutputView' ,imref2d(size(image2)));

diff_image = warpped_im - image2;

[ It_row , It_col ]=size( image1 ); 

mask=zeros(It_row,It_col);        

mask(abs(diff_image) <  thresh) = 0;        
mask(abs(diff_image) >= thresh) = 1;
erode_mask = strel('disk',5);  
dialte_mask = strel('disk',8); 
mask = bwareaopen(mask,2);
mask = imdilate(mask,dialte_mask);
mask = imerode(mask,erode_mask);

        