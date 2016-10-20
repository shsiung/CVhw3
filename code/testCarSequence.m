function testCarSequence()

rect = [60,117,146,152];
data = load('../data/carseq.mat');
data = data.frames;
for i = 1 : 100
   [u, v] = LucasKanadeInverseCompositional(data(:,:,i),data(:,:,i+1),rect);
   rect = round([rect(1)-u, rect(2)-v, rect(3)-u, rect(4)-v]);
end