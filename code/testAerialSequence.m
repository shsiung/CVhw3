function testAerialSequence()

data = load('../data/aerialseq.mat');
data = data.frames;
image_pos = 1;
for i = 1 : size(data,3)-1
    mask = SubtractDominantMotion(data(:,:,i), data(:,:,i+1));  
    if i == 30 || i == 60 || i == 90 || i == 120
       subplot(1,4,image_pos);
       hold on;
       C = imfuse(data(:,:,i),mask,'falsecolor');
       imshow(C);
       image_pos = image_pos + 1;
       str = sprintf('%d',i);
       title(str);
       drawnow;
    end
end