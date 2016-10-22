function testSylvSequence()

rect_basis = [102,62,156,108];
rect_IC = [102,62,156,108];
data = load('../data/sylvseq.mat');
data = data.frames;
bases = load('../data/sylvbases.mat');
bases = bases.bases;
for i = 1 : size(data,3)-1
   clf
   imshow(data(:,:,i));
   hold on;
   rectangle('Position',[rect_basis(1),rect_basis(2),...
                        rect_basis(3)-rect_basis(1),...
                        rect_basis(4)-rect_basis(2)],'EdgeColor','y');
   [u, v] = LucasKanadeBasis(data(:,:,i),data(:,:,i+1),rect_basis,bases);
   
%    hold on;
%    rectangle('Position',[rect_IC(1),rect_IC(2),rect_IC(3)-rect_IC(1),...
%                          rect_IC(4)-rect_IC(2)],'EdgeColor','g');
   %[u1, v1] = LucasKanadeInverseCompositional(data(:,:,i),data(:,:,i+1),rect_IC);
   
   rect_basis = [rect_basis(1)+u, rect_basis(2)+v,...
                 rect_basis(3)+u, rect_basis(4)+v];
%    rect_IC = [rect_IC(1)+u1, rect_IC(2)+v1,... 
%               rect_IC(3)+u1, rect_IC(4)+v1];
   drawnow;
end