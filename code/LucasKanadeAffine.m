% function M = LucasKanadeAffine(It, It1)
% % input - image at time t, image at t+1 
% % output - M affine transformation matrix
% 
% DEBUG = 1;
% 
% %% Precomputation
% % Template gradient - find overlapped region
% It_mask = (double(It)~=0);
% [Gx,Gy] = imgradientxy(double(It));
% 
% steepest_desc = zeros(size(It,1),size(It,2),6);
% 
% [x,y] = meshgrid(1:size(It,2),1:size(It,1));
% steepest_desc(:,:,1) = arrayfun(@(x,y)x,x,y).*Gx;
% steepest_desc(:,:,2) = arrayfun(@(x,y)x,x,y).*Gy;
% steepest_desc(:,:,3) = arrayfun(@(x,y)y,x,y).*Gx;
% steepest_desc(:,:,4) = arrayfun(@(x,y)y,x,y).*Gy;
% steepest_desc(:,:,5) = Gx;
% steepest_desc(:,:,6) = Gy;
% 
% % Hessian 
% H = zeros(6,6);
% for i = 1:6
%     for j = 1:6
%         SD = steepest_desc(:,:,i).*steepest_desc(:,:,j);
%         H(i,j) = sum(SD(:));
%     end
% end
% 
% if (DEBUG)
%     figure(3)
%     imshowpair(uint8(It1),uint8(It),'montage')
% end
% %% Iterate
% p1 = 0;
% p2 = 0;
% p3 = 0; %u
% p4 = 0;
% p5 = 0;
% p6 = 0; %v
% 
% epsilon = 0.15;
% p = [1+p1 p3   p5;...
%      p2   1+p4 p6];
% dP = p;
% 
% if (DEBUG)
%     figure(1)
% end
% M = eye(3);
% SD_curr = steepest_desc;
% while norm(dP) > epsilon
%     % Warp image according to transform matrix
%     tform = affine2d(M');
%     warp_It1 = imwarp(double(It1),tform);
%                     drawnow;
%                     clf
% 
%     % Compute error by substracting the overlap region (by masking)
%     It1_common = zeros(size(It,1),size(It,2));
%     min_width = min(size(It,2),size(warp_It1,2));
%     min_height = min(size(It,1),size(warp_It1,1));
%     It1_common(1:min_height,...
%                1:min_width) = ...
%                warp_It1(1:min_height,...
%                         1:min_width);
%     It1_mask = (It1_common~=0);
% 
%     common_mask = It_mask+It1_mask;
%     It(common_mask~=2) = 0;
%     It1_common(common_mask~=2) = 0; 
%     
%     if (DEBUG)
%          imshowpair(It,It1_common,'montage');
%     end
% 
%     error_image = It1_common-double(It);
%     
%      if (DEBUG)
%                   figure(1)
%                   imshowpair(It1_common,error_image,'montage');
%      end
% 
%     % Per pixel, multiply steepest descent
%     SDparam = zeros(6,1);
%     for i = 1:6
%         SD = steepest_desc(:,:,i);
%         SD(common_mask~=2) = 0;
%         SDlayer = SD.*error_image;
%         SDparam(i) = sum(SDlayer(:));
%     end
% 
%     for i = 1:6
%         SDx = steepest_desc(:,:,i);
%         SDx(common_mask~=2) = 0;
%         for j = 1:6
%             SDy = steepest_desc(:,:,j);            
%             SDy(common_mask~=2) = 0;
%             SD = SDx.*SDy;
%             H(i,j) = sum(SD(:));
%         end
%     end
%     
%     if (DEBUG)
%                figure(4)
%                imshow(H/norm(H))
%     end
% 
%     % Compute delta p
%     dP = pinv()*SDparam;
%     % Update the warp
%     newP = 1/((1+dP(1))*(1+dP(4))-dP(2)*dP(3))*[-dP(1)-dP(1)*dP(4)+dP(2)*dP(3); ...
%                                                 -dP(2);...
%                                                 -dP(3);...
%                                                 -dP(4)-dP(1)*dP(4)+dP(2)*dP(3);...
%                                                 -dP(5)-dP(4)*dP(5)+dP(3)*dP(6);...
%                                                 -dP(6)-dP(1)*dP(6)+dP(2)*dP(5)];
%     
%     p1 = p1+newP(1)+p1*newP(1)+p3*newP(2);
%     p2 = p2+newP(2)+p2*newP(1)+p4*newP(2);
%     p3 = p3+newP(3)+p1*newP(3)+p3*newP(4);
%     p4 = p4+newP(4)+p2*newP(3)+p4*newP(4);
%     p5 = p5+newP(5)+p1*newP(5)+p3*newP(6);
%     p6 = p6+newP(6)+p2*newP(5)+p4*newP(6);
%     
%     p = [1+p1 p3 p5; p2 1+p4 p6];
% %     norm(dP)
%     M = [p ;[0    0    1]];
% end
% 
% if (DEBUG)
%     figure(2)
%     tform = affine2d(inv(M)');
%     warp_It = imwarp(double(It),tform);
%     imshowpair(It1,warp_It,'montage');
% end


function M = LucasKanadeAffine(It, It1)

    It = double(It);
    It1 = double(It1);

    [It_row,It_col]=size(It);%
    x1=1; y1=1;
    x2=It_col; y2=It_row;
    [X,Y]=meshgrid(x1:x2,y1:y2);      
        
    template=It;
    
    [tp_dx ,tp_dy ] = gradient(template);  
   
   
    A=[ X(:).*tp_dx(:)  X(:).*tp_dy(:)  Y(:).*tp_dx(:)  Y(:).*tp_dy(:)  tp_dx(:)   tp_dy(:)] ;
    
    p=zeros(6,1);
    a=0;b=0;c=0;d=0;e=0;f=0;    
    M=[ (1+p(1)) p(3) p(5) ;  p(2) (1+p(4)) p(6) ; 0  0  1]; 
    % ones_tot=ones(size(X(:),1),1)'; 

    step5=A;    
    
    H= step5'*step5;

    invH=pinv(H);

    epsilon=0.05;
    
    error_Total=[];
    %p1=0;p2=0;p3=0;p4=0;p5=0;p6=0;    
    %W=[ (1+p1) p3 p5 ; p2 (1+p4) p6 ];
    Xp=X;
    Yp=Y;
    
    d_p=[10 ;10 ];
    while(sum(abs(d_p))>=epsilon)
        
        %List of points being checked to lie in the template size.
        list_pts_It=[X(:)' ; Y(:)' ; ones(size(X(:),1),1)'];
    
        % Warp Points
        It1_pt_h=M*list_pts_It;
        It1_M_crp_pt = It1_pt_h';
        %It1_M_crp_pt=[ It1_pt_h(1,:)./It1_pt_h(3,:) ; It1_pt_h(2,:)./It1_pt_h(3,:) ;It1_pt_h(3,:)./It1_pt_h(3,:) ];

        %check projection points
        inx=find((It1_M_crp_pt(:,1)<=size(It1,2))& (It1_M_crp_pt(:,1)>=1)&...
                    (It1_M_crp_pt(:,2)<=size(It1,1)) & (It1_M_crp_pt(:,2)>=1) ) ; 
         
        %check for x y replacement !!                
        It1_warpped = interp2(It1,It1_M_crp_pt(inx,1),It1_M_crp_pt(inx,2));         
         
        %ItL = It(:);
        error_image=  -It1_warpped(:) + It(inx);%ItL(inx);        
        
        %only consider values inside the
        %frameab
        step5_curr = step5(inx,:);
        
        d_p=pinv(step5_curr'*step5_curr)*step5_curr'*error_image;
        
        %LK Computation
        %step7=step5_curr'*error_image;
        %delta_p=invH*step7;       
        %W_inv_dp
        
        %{        
        num=1/ ( (1+d_p(1))*(1+d_p(4))- d_p(2)*d_p(3)  );
        mat_2trm=[(-d_p(1) -(d_p(1)*d_p(4)) + (d_p(2)*d_p(3))) ;...
                  -d_p(2);
                  -d_p(3);
                  (-d_p(4) -(d_p(1)*d_p(4)) + d_p(2)*d_p(3) );...
                  (-d_p(5) -(d_p(4)*d_p(5)) + d_p(3)*d_p(6) );...
                  (-d_p(6) -(d_p(1)*d_p(6)) + d_p(2)*d_p(5) ) ];
              
        W_dp_inv=num*mat_2trm;
        
        p=[ (   p(1)+W_dp_inv(1)+p(1)*W_dp_inv(1)+p(3)*W_dp_inv(2) );...
              ( p(2)+W_dp_inv(2)+p(2)*W_dp_inv(1)+p(4)*W_dp_inv(2) );...
              ( p(3)+W_dp_inv(3)+p(1)*W_dp_inv(3)+p(3)*W_dp_inv(4) );...
              ( p(4)+W_dp_inv(4)+p(2)*W_dp_inv(3)+p(4)*W_dp_inv(4) );...
              ( p(5)+W_dp_inv(5)+p(1)*W_dp_inv(5)+p(3)*W_dp_inv(6) );...
              ( p(6)+W_dp_inv(6)+p(2)*W_dp_inv(5)+p(4)*W_dp_inv(6) )...
              ];
        %}
        p = p + d_p;
          M=[1+p(1) p(3) p(5);p(2) p(4)+1 p(6);0 0 1];
        
    end        
