function [u,v] = LucasKanadeInverseCompositional(It, It1, rect)

% input - image at time t, image at t+1, rectangle (top left, bot right coordinates)
% output - movement vector, [u,v] in the x- and y-directions.

%% Precomputation

It = im2double(It);
It1 = im2double(It1);
% Template gradient
[X,Y] = meshgrid(rect(1):rect(3),rect(2):rect(4));
temp_It = interp2(It,X,Y,'linear');
[height,width] = size(temp_It);
[Gx,Gy] = imgradientxy(temp_It);

% Steepest Descent
% Jacobian
% init_dwdp =  [0 0 0 0 1 0; 0 0 0 0 0 1];
steepest_desc = zeros(size(Gx,1),size(Gx,2),6);
steepest_desc(:,:,5) = Gx;
steepest_desc(:,:,6) = Gy;

% Hessian 
H = zeros(6,6);
for i = 1 : size(Gx,2)
    for j = 1: size(Gx,1)
        desc_pixel = reshape(steepest_desc(j,i,:),1,6);
        H = H + desc_pixel' * desc_pixel;
    end
end

%% Iterate
p5 = 0;
p6 = 0; %v

epsilon = 0.01;
p = [1 0 p5; 0 1 p6];
dP = p;
newP = p;

while norm(dP) > epsilon
    % Warp image according to warped window
    % Compute error by subtracting It1 window with It window
    [X1,Y1]=meshgrid(linspace(rect(1)+newP(5),rect(3)+newP(5),width),linspace(rect(2)+newP(6),rect(4)+newP(6),height));
    warp_It1=interp2(It1,X1,Y1,'spline');
    
    error_image = warp_It1 - temp_It;
    
    % Per pixel, multiply steepest descent
    SDparam = zeros(6,1);
    for i = 5:6
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