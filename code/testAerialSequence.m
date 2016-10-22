function testAerialSequence()

data = load('../data/aerialseq.mat');
data = data.frames;
for i = 1 : size(data,3)-1
   clf
   mask = SubtractDominantMotion(data(:,:,i), data(:,:,i+1));  
   imshow(mask);
end