function testUSSeqAffine()

data = load('../data/usseq.mat');
data = data.frames;
rects = load('../results/usseqrects.mat');
rects = rects.usseqrects;
image_pos = 1;
rect = rects(1,:);
for i = 1 : size(data,3)-1
    width = rect(3)-rect(1);
    height = rect(4)-rect(2);
    A = (0:width) + rect(1);
    B = (0:height) + rect(2);
    [X,Y] = meshgrid(A, B);
    It = interp2(double(data(:,:,i)),X,Y,'spline');
    It1 = interp2(double(data(:,:,i+1)),X,Y,'spline');

    if i == 4 || i == 24 || i == 49 || i == 74 || i == 99
       subplot(1,5,image_pos);
       hold on;
       C = imfuse(data(:,:,1),mask,'falsecolor');
       hold on;
       imshow(C)
       image_pos = image_pos + 1;
       str = sprintf('%d',i+1);
       title(str);
       drawnow;
    end
    mask = SubtractDominantMotion(It, It1); 
    rect = rects(i+1,:);
end