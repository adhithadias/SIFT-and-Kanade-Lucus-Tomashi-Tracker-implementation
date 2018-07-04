close all

imagefiles = dir('graf/*.ppm');      
nfiles = length(imagefiles);    % Number of files found
number_of_points = 6;

for ii=1:nfiles
  disp( ['Reading image ' num2str(ii)] );
  currentimage = imread(['graf/' imagefiles(ii).name]);
  #images{ii} = rgb2gray(currentimage);
  #images{ii} = im2single(images{ii});
  images{ii} = currentimage;
  #figure, imshow(images{ii})
end

subplot(1,2,1)
imagesc(images{1}), hold on, axis image
[x1,y1] = ginput(number_of_points);
plot(x1, y1, 'g'); drawnow; hold off;

subplot(1,2,2)
imagesc(images{2}), hold on, axis image
[x2,y2] = ginput(number_of_points);
plot(x2, y2, 'g'); drawnow; hold off;

a2 = [x2'; y2'];
a1 = [x1'; y1'; ones(1, numel(x1))]; 

#a2 = trans_mat * a1;
trans_mat = a2 * pinv(a1);
trans_mat(3,:) = [0, 0, 1];

a2
trans_mat * a1

#J = imwarp(images{1},trans_mat);
J = imperspectivewarp(images{1}, trans_mat, "bicubic", "same");
figure,
imagesc(J), hold on, axis image