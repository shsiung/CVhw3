function testCarSequence()

rect = [60,117,146,152];
data = load('../data/carseq.mat');
data = data.frames;
carseqrects  = zeros(size(data,3),4);
carseqrects(1,:) = rect;
image_pos = 1;
for i = 1 : size(data,3)-1
   if i == 2 || i == 100 || i == 200 || i == 300 || i == 400
       subplot(1,5,image_pos);
       hold on;
       imshow(data(:,:,i));
       drawnow;
       rectangle('Position',[rect(1),rect(2),rect(3)-rect(1),rect(4)-rect(2)],'EdgeColor','g','LineWidth',2);
       image_pos = image_pos + 1;
       str = sprintf('%d (%f milliseconds)',i, elapsedTime);
       title(str);
   end
   tic
   [u, v] = LucasKanadeInverseCompositional(data(:,:,i),data(:,:,i+1),rect);
   rect = [rect(1)+u, rect(2)+v, rect(3)+u, rect(4)+v];
   carseqrects(i+1,:) = rect;
   elapsedTime = toc * 1000;
end

save('carseqrects.mat','carseqrects');