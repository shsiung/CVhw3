function testUltrasoundSequence()

rect = [255,105,310,170];
data = load('../data/usseq.mat');
data = data.frames;
for i = 1 : size(data,3)-1
    i
   imshow(data(:,:,i));
   hold on;
   rectangle('Position',[rect(1),rect(2),rect(3)-rect(1),rect(4)-rect(2)]);
   drawnow;
   [u, v] = LucasKanadeInverseCompositional(data(:,:,i),data(:,:,i+1),rect);
   
   rect = [rect(1)-u, rect(2)-v, rect(3)-u, rect(4)-v];
end