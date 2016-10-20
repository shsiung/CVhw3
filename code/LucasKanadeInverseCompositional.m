function [u,v] = LucasKanadeInverseCompositional(It, It1, rect)

% input - image at time t, image at t+1, rectangle (top left, bot right coordinates)
% output - movement vector, [u,v] in the x- and y-directions.

It = im2double(It);
It1 = im2double(It1);
%% Precomputation

% Template gradient
[Gx, Gy] = imgradientxy(It,'intermediate');
whole_rect = round(rect);
Gx = imtranslate(Gx,[whole_rect(1)-rect(1),whole_rect(2)-rect(2)]);
Gy = imtranslate(Gy,[whole_rect(1)-rect(1),whole_rect(2)-rect(2)]);
Gx = Gx(whole_rect(2):whole_rect(4), whole_rect(1):whole_rect(3));
Gy = Gy(whole_rect(2):whole_rect(4), whole_rect(1):whole_rect(3));

% Steepest Descent
% Jacobian
% init_dwdp =  [0 0 0 0 1 0; 0 0 0 0 0 1];
steepest_desc = zeros(whole_rect(4)-whole_rect(2)+1,whole_rect(3)-whole_rect(1)+1,6);
steepest_desc(:,:,5) = Gx;
steepest_desc(:,:,6) = Gy;

% Hessian 
H = zeros(6,6);
for i = 1 : whole_rect(3)-whole_rect(1)+1
    for j = 1: whole_rect(4)-whole_rect(2)+1
        desc_pixel = reshape(steepest_desc(j,i,:),1,6);
        H = H + desc_pixel' * desc_pixel;
    end
end

%% Iterate
p5 = 0;
p6 = 0; %v

epsilon = 0.01;
p = [1 0 p5; 0 1 p6];
temp_It = imtranslate(It,[whole_rect(1)-rect(1),whole_rect(2)-rect(2)]);
temp_It = temp_It(whole_rect(2):whole_rect(4), whole_rect(1):whole_rect(3));
dP = p;
newP = p;

while norm(dP) > epsilon
    % Warp image according to warped window
    warp_It1 = imtranslate(It1,[-newP(5),-newP(6)]);
    warp_It1 = warp_It1(whole_rect(2):whole_rect(4), whole_rect(1):whole_rect(3));
    % Compute error by subtracting It1 window with It window
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
    p5 = p5+dP(5);
    p6 = p6+dP(6);
    newP = [0;0;0;0;-p5;-p6];
    u = newP(5);
    v = newP(6);
    
end


