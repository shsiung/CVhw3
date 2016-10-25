function M = LucasKanadeAffine(It, It1)

% input - image at time t, image at t+1
% output - transformation matrix M

It = double(It);
It1 = double(It1);

[It_row,It_col]=size(It);
[X,Y]=meshgrid(1:It_col,1:It_row);      

[Gx ,Gy ] = gradient(It);  
SD=[ X(:).*Gx(:)  X(:).*Gy(:)  Y(:).*Gx(:)  Y(:).*Gy(:)  Gx(:)   Gy(:)] ;

p=zeros(6,1);
M=[ (1+p(1)) p(3) p(5) ;  p(2) (1+p(4)) p(6) ; 0  0  1]; 

epsilon=0.05;

dP=[10 ;10 ];
while(sum(abs(dP))>=epsilon)

    % List of points being checked to lie in the template size.
    It_pts=[X(:)' ; Y(:)' ; ones(size(X(:),1),1)'];

    % Warp Points
    warped_It_pts=M*It_pts;
    crop_It_pts = warped_It_pts';

    % Check common points
    common_idx=find((crop_It_pts(:,1)<=size(It1,2))& (crop_It_pts(:,1)>=1)&...
                    (crop_It_pts(:,2)<=size(It1,1)) & (crop_It_pts(:,2)>=1) ) ; 

    % Check for x y replacement !!                
    It1_warpped = interp2(It1,crop_It_pts(common_idx,1),crop_It_pts(common_idx,2));         

    error_image = -It1_warpped(:) + It(common_idx);       

    % Only consider values inside the
    SD_curr = SD(common_idx,:);

    % Calculate steepest descent as well as the hession here.
    dP=pinv(SD_curr'*SD_curr)*SD_curr'*error_image;

    p = p + dP;
    M=[1+p(1) p(3) p(5);p(2) p(4)+1 p(6);0 0 1];

end  
