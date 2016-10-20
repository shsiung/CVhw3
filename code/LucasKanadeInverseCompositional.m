function [u,v] = LucasKanadeInverseCompositional(It, It1, rect)

% input - image at time t, image at t+1, rectangle (top left, bot right coordinates)
% output - movement vector, [u,v] in the x- and y-directions.
It = im2double(It);
It1 = im2double(It1);
%% Precomputation
% Jacobian
init_dwdp = [0 0 1 0 0 0; 0 0 0 0 0 1];
% Template gradient
[Gx, Gy] = imgradientxy(It,'intermediate');

H = zeros(6,6);
steepest_desc = zeros(rect(4)-rect(2)+1,rect(3)-rect(1)+1,6);
for i = 1 : rect(3)-rect(1)+1
    for j = 1: rect(4)-rect(2)+1
        % Steepest Descent
        steepest_desc(j,i,:) = [Gx(j,i), Gy(j,i)] * init_dwdp;
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

epsilon = 0.08;
p = [1+p1 p2 p3; p4 1+p5 p6];
W = @(x,p) p*[x,1]';
temp_It = It(rect(2):rect(4), rect(1):rect(3));
dP = p;
newP = p;
rec_corner = [rect(1) rect(2);...
              rect(1) rect(4);...
              rect(3) rect(2);...
              rect(3) rect(4)];
newRec = [rec_corner(1,1), rec_corner(1,2), rec_corner(4,1), rec_corner(4,2)];
[M, N]=size(It1); %Assuming square matrix.
[xx, yy]=meshgrid(1:N,1:M); %xx,yy are both outputs of meshgrid (so called plaid matrices).

%=================================== Plotting for debugging ==============
subplot(2,2,1)
imshow(temp_It);
subplot(2,2,2)
orig_It1=It1(rect(2):rect(4), rect(1):rect(3));
imshow(orig_It1);
%=================================== Plotting for debugging ==============
% error_image = [2 2; 2 2];
while norm(dP) > epsilon
    % Warp image according to warped window
    warp_It1=imtranslate(It1,[-newP(1,3),-newP(2,3)]);
%     warp_It1=interp2(xx,yy,It1,xx+newP(1,3),yy+newP(2,3));
    warp_It1=warp_It1(rect(2):rect(4), rect(1):rect(3));
    % Compute error by subtracting It1 window with I;t window
    error_image = warp_It1 - temp_It;
    
    %imshowpair(temp_It,error_image,'montage');
    % Per pixel, multiply steepest descent
    SDparam = zeros(6,1);

    for i = 1: rect(3)-rect(1)+1
        for j = 1 : rect(4)-rect(2)+1
            % Steepest Descent
            SDparam = SDparam + reshape(steepest_desc(j,i,:),6,1)*error_image(j,i);
        end
    end

    % Compute delta p
    dP = pinv(H)*SDparam;
    % Update the warp
    p1 = p1+dP(1)+p1*dP(1)+p3*dP(2);
    p2 = p2+dP(2)+p2*dP(1)+p4*dP(2);
    p3 = p3+dP(3)+p1*dP(3)+p3*dP(4);
    p4 = p4+dP(4)+p2*dP(3)+p4*dP(4);
    p5 = p5+dP(5)+p1*dP(5)+p3*dP(6);
    p6 = p6+dP(6)+p2*dP(5)+p4*dP(6);
    newP = 1/((1+p1)*(1+p4)-p2*p3)*[-p1-p1*p4+p2*p3, ...
                                           -p2,...
                                           -p3;...
                                         -p4-p1*p4+p2*p3,...
                                         -p5-p4*p5+p3*p6,...
                                         -p6-p1*p6+p2*p5];
    subplot(2,2,4);
    imshow(warp_It1);
    subplot(2,2,3);
    imshow(error_image);
    drawnow;
    %rec_corner = arrayfun(@(x,p) W(x,newp), rec_corner, newP);
    %newRec = [rec_corner(1,1), rec_corner(1,2), rec_corner(4,1), rec_corner(4,2)];
    [sum(sum(error_image)) norm(dP)]
end
u = p3;
v = p6;
