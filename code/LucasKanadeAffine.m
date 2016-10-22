function M = LucasKanadeAffine(It, It1)
% input - image at time t, image at t+1 
% output - M affine transformation matrix

It = im2double(It);
It1 = im2double(It1);
%% Precomputation
% Template gradient - find overlapped region
It_mask = (It~=0);
[Gx,Gy] = imgradientxy(It);

init_dwdp =@(x,y)  [x/size(It1,2) 0 y/size(It1,1) 0 1 0;...
                    0 x/size(It1,2) 0 y/size(It1,1) 0 1]; % Jacobian
steepest_desc = zeros(size(It,2),size(It,1),6);

H = zeros(6,6);
for i = 1 : size(It,2)
    for j = 1: size(It,1)
        % Steepest Descent
        steepest_desc(:,:,i) = [Gx(j,i), Gy(j,i)].*init_dwdp(j,i);
        desc_pixel = reshape(steepest_desc(j,i,:),1,6);
       
         % Hessian 
        H = H + desc_pixel' * desc_pixel;
    end
end

%% Iterate
p1 = 0;
p2 = 0;
p3 = 0; %u
p4 = 0;
p5 = 0;
p6 = 0; %v

epsilon = 0.01;
p = [1+p1 p3 p5; p2 1+p4 p6];
dP = p;
while norm(dP) > epsilon
    % Warp image according to transform matrix
    tform = affine2d([p;[0,0,1]]');
    warp_It1 = imwarp(It1,tform);
  
    % Compute error by substracting the overlap region (by masking)
    It1_common = zeros(size(It,2),size(It,1));
    It1_common(1:min(size(It,2),size(warp_It1,2)),...
               1:min(size(It,1),size(warp_It1,1))) = ...
               warp_It1(1:min(size(It,2),size(warp_It1,2)),...
                        1:min(size(It,1),size(warp_It1,1)));
    
    common_mask = It_mask+It_common;
    It(common_mask~=2) = 0;
    It1_common(common_mask~=2) = 0; 
    error_image = It1_common-It;
    
    % Per pixel, multiply steepest descent
    SDparam = zeros(6,1);
    
    for i = 1:6
        SD = steepest_desc(:,:,1);
        SD(common_mask~=2) = 0;
        SDlayer = SD.*error_image;
        SDparam(i) = sum(SDlayer(:));
    end

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

end
M = [1+p1 p3   p5;...
     p2   1+p4 p6;...
     0    0    1];
