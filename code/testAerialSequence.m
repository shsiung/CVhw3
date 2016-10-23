function testAerialSequence()

data = load('../data/aerialseq.mat');
data = data.frames;
for i = 1 : size(data,3)-1
   i
   clf;
   mask = SubtractDominantMotion(data(:,:,i), data(:,:,i+1));  
   C = imfuse(mask,data(:,:,i),'blend');
   imshow(C);
   drawnow;
end