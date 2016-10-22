function M = LucasKanadeAffine(It, It1)
% input - image at time t, image at t+1 
% output - M affine transformation matrix

%% Precomputation
% Template gradient - find overlapped region
It_mask = (double(It)~=0);
[Gx,Gy] = imgradientxy(double(It));

steepest_desc = zeros(size(It,1),size(It,2),6);
[x,y] = meshgrid(1:size(It,2),1:size(It,1));
steepest_desc(:,:,1) = arrayfun(@(x,y)x,x,y).*Gx;
steepest_desc(:,:,2) = arrayfun(@(x,y)x,x,y).*Gy;
steepest_desc(:,:,3) = arrayfun(@(x,y)y,x,y).*Gx;
steepest_desc(:,:,4) = arrayfun(@(x,y)y,x,y).*Gy;
steepest_desc(:,:,5) = Gx;
steepest_desc(:,:,6) = Gy;
%imshowpair(steepest_desc(:,:,3),steepest_desc(:,:,4),'montage');
%imshow(steepest_desc(:,:,5))
% figure
% imshow(steepest_desc(:,:,1))
% figure
% imshow(steepest_desc(:,:,2))
% figure
% imshow(steepest_desc(:,:,3))
% figure
% imshow(steepest_desc(:,:,4))
% figure
% imshow(steepest_desc(:,:,5))
% figure
% imshow(steepest_desc(:,:,6))
% figure

% Hessian 
H = zeros(6,6);
for i = 1:6
    for j = 1:6
        SD = steepest_desc(:,:,i)'*steepest_desc(:,:,j);
        H(i,j) = sum(SD(:));
    end
end
figure(3)
imshowpair(uint8(It1),uint8(It),'montage')
figure(4)
imshow(uint8(H));
%% Iterate
p1 = 0;
p2 = 0;
p3 = 0; %u
p4 = 0;
p5 = 0;
p6 = 0; %v

epsilon = 0.01;
p = [1+p1 p3   p5;...
     p2   1+p4 p6];
dP = p;
figure(1)
iter = 0;
M = eye(3);
while norm(dP) > epsilon
    iter = iter + 1
    % Warp image according to transform matrix
    tform = affine2d(M');
    warp_It1 = imwarp(double(It1),tform);
    hold on;
    drawnow;
        clf

    % Compute error by substracting the overlap region (by masking)
    It1_common = zeros(size(It,1),size(It,2));
    min_width = min(size(It,2),size(warp_It1,2));
    min_height = min(size(It,1),size(warp_It1,1));
    It1_common(1:min_height,...
               1:min_width) = ...
               warp_It1(1:min_height,...
                        1:min_width);
    It1_mask = (It1_common~=0);

    common_mask = It_mask+It1_mask;
    It(common_mask~=2) = 0;
    It1_common(common_mask~=2) = 0; 
         %imshowpair(It,It1_common,'montage')

    error_image = It1_common-double(It);
                 figure(1)
                 imshowpair(It1_common,error_image,'montage');

    % Per pixel, multiply steepest descent
    SDparam = zeros(6,1);
    
            %figure
    for i = 1:6
        SD = steepest_desc(:,:,i);
        SD(common_mask~=2) = 0;
            %subplot(2,3,i)
            %imshow(SD)
        SDlayer = SD'*error_image;
        SDparam(i) = sum(SDlayer(:));
    end
    
    for i = 1:6
        for j = 1:6
            SD = steepest_desc(:,:,i);
            SD(common_mask~=2) = 0;
            SD = SD'*SD;
            H(i,j) = sum(SD(:));
        end
    end
    
    figure(4)
    imshow(H)

    % Compute delta p
    dP = pinv(H)*SDparam;
    
    % Update the warp
    newP = 1/((1+dP(1))*(1+dP(4))-dP(2)*dP(3))*[-dP(1)-dP(1)*dP(4)+dP(2)*dP(3); ...
                                                -dP(2);...
                                                -dP(3);...
                                                -dP(4)-dP(1)*dP(4)+dP(2)*dP(3);...
                                                -dP(5)-dP(4)*dP(5)+dP(3)*dP(6);...
                                                -dP(6)-dP(1)*dP(6)+dP(2)*dP(5)];
    
    p1 = p1+newP(1)+p1*newP(1)+p3*newP(2);
    p2 = p2+newP(2)+p2*newP(1)+p4*newP(2);
    p3 = p3+newP(3)+p1*newP(3)+p3*newP(4);
    p4 = p4+newP(4)+p2*newP(3)+p4*newP(4);
    p5 = p5+newP(5)+p1*newP(5)+p3*newP(6);
    p6 = p6+newP(6)+p2*newP(5)+p4*newP(6);
    
    p = [1+p1 p3 p5; p2 1+p4 p6];
    norm(dP)
    M = [p ;[0    0    1]];
end


figure(2)
tform = affine2d(M');
warp_It1 = imwarp(double(It1),tform);
imshowpair(It,warp_It1,'montage');
