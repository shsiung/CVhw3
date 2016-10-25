function testUltrasoundSequence()

rect = [255; 105; 310; 170];
data = load('../data/usseq.mat');
data = data.frames;
usseqrects  = zeros(size(data,3),4);
usseqrects(1,:) = rect;
image_pos = 1;
for i = 1 : size(data,3)-1
   if i == 4 || i == 24 || i == 49 || i == 74 || i == 99
       subplot(1,5,image_pos);
       hold on;
       imshow(data(:,:,i));
       drawnow;
       rectangle('Position',[rect(1),rect(2),rect(3)-rect(1),rect(4)-rect(2)],'EdgeColor','g','LineWidth',2);
       image_pos = image_pos + 1;
       str = sprintf('%d (%f milliseconds)',i+1, elapsedTime);
       title(str);
   end
   tic
   [u, v] = LucasKanadeInverseCompositional(data(:,:,i),data(:,:,i+1),rect);
   rect = [rect(1)+u, rect(2)+v, rect(3)+u, rect(4)+v];
   usseqrects(i+1,:) = rect;
   elapsedTime = toc * 1000;
end

save('usseqrects.mat','usseqrects');