function M = LucasKanadeAffine(It, It1)
% input - image at time t, image at t+1 
% output - M affine transformation matrix

It = im2double(It);
It1 = im2double(It1);
%% Precomputation

% Template gradient - find overlapped region

temp_It = It;
It_mask = (temp_It~=0);

[Gx,Gy] = imgradientxy(temp_It);

init_dwdp =@(x,y)  [x 0 y 0 1 0; 0 x 0 y 0 1]; % Jacobian
steepest_desc = zeros(whole_rect(4)-whole_rect(2)+1,whole_rect(3)-whole_rect(1)+1,6);

H = zeros(6,6);
for i = 1 : whole_rect(3)-whole_rect(1)+1
    for j = 1: whole_rect(4)-whole_rect(2)+1
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
newP = p;

while norm(dP) > epsilon
    % Warp image according to transform matrix
    tform = affine2d([p;[0,0,1]]');
    warp_It1 = imwarp(It1,tform);
    It1_mask = (warp_It1~=0);
    % Compute error by substracting the overlap region
    max_width = max(size(It1_mask,1),size(It_mask,1));
    max_height = max(size(It1_mask,2),size(It_mask,2));
    It1_common = zeros(max_height,max_width);
    It1_common(1:size(It1_mask,2),1:size(It1_mask,1)) = It1_mask;
    It_common = zeros(max_height,max_width);
    It_common(1:size(It_mask,2),1:size(It_mask,1)) = It_mask;
    
    error_image = It1_common+It_common;
    error_image = (error_image==2);
    
    error_image = warp_It1 - temp_It;
    
    % Per pixel, multiply steepest descent
    SDparam = zeros(6,1);

    for i = 1:6
        SDlayer = steepest_desc(:,:,i).*error_image;
        SDparam(i) = sum(SDlayer(:));
    end

    % Compute delta p
    dP = pinv(H)*SDparam;
    
    % Update the warp
    newP = 1/((1+dP(1))*(1+dP(4))-p2*p3)*[-dP1-dP1*p4+p2*p3; ...
                                    -p2;...
                                    -p3;...
                                    -p4-p1*p4+p2*p3;...
                                    -p5-p4*p5+p3*p6;...
                                    -p6-p1*p6+p2*p5];
    
    p1 = p1+newP(1)+p1*newP(1)+p3*newP(2);
    p2 = p2+newP(2)+p2*newP(1)+p4*newP(2);
    p3 = p3+newP(3)+p1*newP(3)+p3*newP(4);
    p4 = p4+newP(4)+p2*newP(3)+p4*newP(4);
    p5 = p5+newP(5)+p1*newP(5)+p3*newP(6);
    p6 = p6+newP(6)+p2*newP(5)+p4*newP(6);
    
end

M = [1+newP(1) newP(3) newP(5);...
     newP(2) 1+newP(4) newP(6);...
     0       0         1     ];
