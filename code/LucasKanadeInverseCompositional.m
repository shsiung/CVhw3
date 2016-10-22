% function [u,v] = LucasKanadeInverseCompositional(It, It1, rect)
% % input - image at time t, image at t+1, rectangle (top left, bot right coordinates)
% % output - movement vector, [u,v] in the x- and y-directions.
% 
% %% Precomputation
% width = rect(3)-rect(1);
% height = rect(4)-rect(2);
% 
% % Template gradient
% A = (0:width) + rect(1);
% B = (0:height) + rect(2);
% [X,Y] = meshgrid(A, B);
% temp_It = interp2(double(It),X,Y,'spline');
% [Gx,Gy] = imgradientxy(temp_It);
% 
% % Steepest Descent
% % Jacobian
% % init_dwdp =  [0 0 0 0 1 0; 0 0 0 0 0 1];
% steepest_desc = zeros(size(Gx,1),size(Gx,2),6);
% steepest_desc(:,:,5) = Gx;
% steepest_desc(:,:,6) = Gy;
% 
% % Hessian 
% H = zeros(6,6);
% for i = 5:6
%     for j = 5:6
%         SD = steepest_desc(:,:,i)'*steepest_desc(:,:,j);
%         H(i,j) = sum(SD(:));
%     end
% end
% 
% %% Iterate
% p5 = 0;
% p6 = 0; %v
% 
% epsilon = 0.01;
% dP = [1 0 p5; 0 1 p6];
% 
% while norm(dP) > epsilon
%     
%     % Update the warp
%     p5 = p5-dP(5);
%     p6 = p6-dP(6); 
%     
%     % Warp image according to warped window
%     % Compute error by subtracting It1 window with It window
%     newA = A+p5;
%     newB = B+p6;
%     [X1,Y1] = meshgrid(newA, newB);
%     warp_It1 = interp2(double(It1),X1,Y1,'spline');
% 
%     error_image = warp_It1 - temp_It;
%     
%     % Multiply steepest descent
%     SDparam = zeros(6,1);
%     for i = 5:6
%         SDlayer = steepest_desc(:,:,i)'*error_image;
%         SDparam(i) = sum(SDlayer(:));
%     end
% 
%     % Compute delta p
%     dP = pinv(H)*SDparam;
% end
% 
% u = p5;
% v = p6;

function [u,v] = LucasKanadeInverseCompositional(It, It1, rect)

% input - image at time t, image at t+1, rectangle (top left, bot right coordinates)
% output - movement vector, [u,v] in the x- and y-directions.

%It = im2double(It);
%It1 = im2double(It1);



tol = 0.1;
max_iter = 100;

rect_new1 = [rect(1) rect(2)]';
err = tol;
iter = 0;
rect_w = rect(3) - rect(1);
rect_h = rect(4) - rect(2);
A = 0:rect_w ;
B =  0:rect_h;

[XI,YI] = meshgrid( A + rect(1), B+ rect(2) );
T = interp2(double(It),XI,YI,'spline');   


N = numel(T);
X = 1:size(T,1);
Y = 1:size(T,2);


[Tx,Ty] = gradient(double(T));
Tj = [reshape(Tx,N,1) reshape(Ty,N,1)];
p = zeros(2,1);
Wj = eye(2);

TjWj = zeros(1,2,N);
H = zeros(2,2);
for i = 1:N
    TjWj(:,:,i) = Tj(i,:) * Wj;
    H = H + TjWj(:,:,i)'*TjWj(:,:,i);
end
H_inv = pinv(H);

while err >= tol && iter < max_iter
    iter = iter+1;

    rect_new1 = [rect(1) rect(2)]' + p;% M*rect_new1;
    rect_new2 = rect_new1 + [rect_w rect_h]'; %M*rect_new2;
    if(rect_new1(1) < 1 || rect_new2(1)<1 || rect_new1(2) > size(It,2) || rect_new2(2) > size(It,2))
        break;
    end
    
    [XI,YI] = meshgrid( A + rect_new1(1), B+ rect_new1(2) );
    newI = interp2(double(It1),XI,YI,'spline');

    
%      subplot(1,2,1);
%     imshow(uint8(T));
%     subplot(1,2,2);
%     imshow(uint8(newI));
    
    Idiff = newI - T;
    s = 0;
    
   % err = sum(sum((Idiff).^2));
    
%     subplot(1,2,1);
%     imshow(uint8(T));
%     subplot(1,2,2);
%     imshow(uint8(newI));

    for i = 1:N
        s = s+ TjWj(:,:,i)'*Idiff(i);
    end
    
    pd = H_inv*s; %delta p
    pd_inv = -1*pd;
    p = p + pd_inv;

  % norm(pd)
    if norm(pd) <0.00001
         break;
    end

end

u = p(1);
v = p(2);






