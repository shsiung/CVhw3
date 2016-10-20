function testCarSequence()

rect = [60,117,146,152];
data = load('../data/carseq.mat');
data = data.frames;
for i = 1 : size(data,3)-1
   
   imshow(data(:,:,i));
   hold on;
   rectangle('Position',[rect(1),rect(2),rect(3)-rect(1),rect(4)-rect(2)]);
   drawnow;
   [u, v] = LucasKanadeInverseCompositional(data(:,:,i),data(:,:,i+1),rect);
   rect = [rect(1)+u, rect(2)+v, rect(3)+u, rect(4)+v];
end