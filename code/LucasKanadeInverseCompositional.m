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

% function [u,v] = LucasKanadeInverseCompositional3(It, It1, rect)
% % It is the template
% % input - image at time t, image at t+1, rectangle (top left, bot right coordinates)
% % output - movement vector, [u,v] in the x- and y-directions.
% 
% 
% It=im2double(It);It1=im2double(It1);
% 
% [X,Y]=meshgrid(rect(1):rect(3),rect(2):rect(4));
% Itwin=interp2(It,X,Y,'linear');
% [height,width]=size(Itwin);
% %figure();imshow(Itwin);
% %[gradTx,gradTy]=gradient(double(Itwin));
% [gradTx,gradTy]=imgradientxy(Itwin);
% Hessian=zeros(6);
% Hessian(5,5)=sum(sum(gradTx.'*gradTx));
% Hessian(5,6)=sum(sum(gradTx.'*gradTy));
% Hessian(6,5)=sum(sum(gradTy.'*gradTx));
% Hessian(6,6)=sum(sum(gradTy.'*gradTy));
% invHessian=pinv(Hessian);
% % trans_Hessian=[sum(sum(gradTx'*gradTx)),sum(sum(gradTx'*gradTy));sum(sum(gradTy'*gradTx)),sum(sum(gradTy'*gradTy))];
% % invtrans_Hessian=inv(trans_Hessian);
% 
% u=0;v=0; %Initial value of u,v
% 
% %figure();
% while (true)
%     
%     % warping and interpolate
% %     It1_warp=imtranslate(It1,[u,v],'linear'); % move image to fit window  
% %     It1win=It1_warp(rect(2):rect(4),rect(1):rect(3));
%     [X1,Y1]=meshgrid(linspace(rect(1)-u,rect(3)-u,width),linspace(rect(2)-v,rect(4)-v,height));
%     It1win=interp2(It1,X1,Y1,'spline');
%     %imshow(It1win);
%     
%     error_mat=It1win-Itwin;
%     %error=sum(sum(abs(error_mat)));
%     
%     
%     intermediate_sol=zeros(6,1);
%     intermediate_sol(5,1)=sum(sum(gradTx.'*error_mat));
%     intermediate_sol(6,1)=sum(sum(gradTy.'*error_mat));
%    
%     delta_p=invHessian*intermediate_sol;    
%     u=u+delta_p(5,1);v=v+delta_p(6,1);
% 
% %     trans_b=[sum(sum(gradTx'*error_mat));sum(sum(gradTy'*error_mat))];
% %     delta_p=trans_Hessian\trans_b;
% %     u=u+delta_p(1);v=v+delta_p(2);
%     
%     % 
%     if (norm(delta_p)<0.01); %disp([u,v,norm(delta_p)]);
%         break;
%     else  %disp([u,v,norm(delta_p)]);
%     end;
%       
% end