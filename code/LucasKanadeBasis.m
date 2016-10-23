function [u,v] = LucasKanadeBasis(It, It1, rect, bases)

% input - image at time t, image at t+1, rectangle (top left, bot right
% coordinates), bases 
% output - movement vector, [u,v] in the x- and y-directions.

%% Precomputation
width = size(bases,2)-1;
height = size(bases,1)-1;

% Template gradient
A = (0:width) + rect(1);
B = (0:height) + rect(2);
[X,Y] = meshgrid(A, B);
temp_It = interp2(double(It),X,Y,'spline');
[Gx,Gy] = imgradientxy(temp_It);

% Steepest Descent
% Jacobian
% init_dwdp =  [0 0 0 0 1 0; 0 0 0 0 0 1];
steepest_desc = zeros(size(Gx,1),size(Gx,2),6);
steepest_desc(:,:,5) = Gx;
steepest_desc(:,:,6) = Gy;

% Weighted bases
for m = 1 : size(bases,3)
    weight = zeros(1,1,6);
   for i = 5:6
    weight(:,:,i) = sum(sum(steepest_desc(:,:,i).*bases(:,:,m)));
   end
   for i = 1 : 6
    steepest_desc(:,:,i) = steepest_desc(:,:,i) - weight(:,:,i) * bases(:,:,m);
   end
end
% Hessian 
H = zeros(6,6);
for i = 5:6
    for j = 5:6
        SD = steepest_desc(:,:,i).*steepest_desc(:,:,j);
        H(i,j) = sum(SD(:));
    end
end

%% Iterate
p5 = 0;
p6 = 0; %v

epsilon = 0.01;
p = [1 0 p5; 0 1 p6];
dP = p;
lambda = zeros(size(bases,3),1);

while norm(dP) > epsilon
    % Warp image according to warped window
    % Compute error by subtracting It1 window with It window
    newA = A+p5;
    newB = B+p6;
    [X1,Y1] = meshgrid(newA, newB);
    warp_It1 = interp2(double(It1),X1,Y1,'spline');

    error_image = warp_It1 - temp_It - ...
        sum(bsxfun(@times,bases,reshape(lambda,[1 1 size(bases,3)])),3);
    
    % Per pixel, multiply steepest descent
    SDparam = zeros(6,1);
    
    for i = 5:6
        SDlayer = steepest_desc(:,:,i).*error_image;
        SDparam(i) = sum(SDlayer(:));
    end

    % Compute delta p
    dP = pinv(H)*SDparam;
    
    % Update the warp
    p5 = p5-dP(5);
    p6 = p6-dP(6);
    
    norm(dP);
    
    % Update lambda
    for i = 1 : size(bases,3)
        lambda(i) = sum(sum(bases(:,:,i) .* (warp_It1 - temp_It)));
    end
end

u = p5;
v = p6;