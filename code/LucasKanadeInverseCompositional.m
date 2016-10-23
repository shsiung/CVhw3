function [u,v] = LucasKanadeInverseCompositional(It, It1, rect)
% input - image at time t, image at t+1, rectangle (top left, bot right coordinates)
% output - movement vector, [u,v] in the x- and y-directions.

%% Precomputation
width = rect(3)-rect(1);
height = rect(4)-rect(2);

% Template gradient
A = (0:width) + rect(1);
B = (0:height) + rect(2);
[X,Y] = meshgrid(A, B);
temp_It = interp2(double(It),X,Y,'linear');
[Gx,Gy] = imgradientxy(temp_It,'central');

% Steepest Descent
% Jacobian
steepest_desc = zeros(size(Gx,1),size(Gx,2),6);
steepest_desc(:,:,5) = Gx;
steepest_desc(:,:,6) = Gy;

% Hessian 
H = zeros(6,6);
for i = 5:6
    for j = 5:6
        SD = steepest_desc(:,:,i).*steepest_desc(:,:,j);
        H(i,j) = sum(SD(:));
    end
end

invH = pinv(H);

%% Iterate
p5 = 0;
p6 = 0; %v

epsilon = 0.2;
dP = [1 0 p5; 0 1 p6];

while norm(dP) > epsilon
    
    % Update the warp
    p5 = p5-dP(5);
    p6 = p6-dP(6); 
    
    % Warp image according to warped window
    % Compute error by subtracting It1 window with It window
    newA = A+p5;
    newB = B+p6;
    [X1,Y1] = meshgrid(newA, newB);
    warp_It1 = interp2(double(It1),X1,Y1,'linear');

    error_image = warp_It1 - temp_It;
    
    % Multiply steepest descent
    SDparam = zeros(6,1);
    for i = 5:6
        SDlayer = steepest_desc(:,:,i).*error_image;
        SDparam(i) = sum(SDlayer(:));
    end
    
    % Compute delta p
    dP = invH*SDparam;
end

u = p5;
v = p6;