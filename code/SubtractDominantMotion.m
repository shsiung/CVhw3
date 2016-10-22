function mask = SubtractDominantMotion(image1, image2)

% input - image1 and image2 form the input image pair
% output - mask is a binary image of the same size

M = LucasKanadeAffine(image1,image2);
thresh = 80;
tform = affine2d(M');
warp_im1 = imwarp(double(image1),tform);

result_im1 = warp_im1(1:size(image1,1),1:size(image1,2));

error = abs(double(image2)-result_im1);

error(error<=thresh)=0;
error(error>thresh)=1;

mask = error;


imshow(mask)