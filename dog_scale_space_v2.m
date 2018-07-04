clear all
close all

imagefiles = dir('graf/*.ppm');      
nfiles = length(imagefiles);    % Number of files found

for ii=1:nfiles
  disp(ii)
  currentimage = imread(['graf/' imagefiles(ii).name]);
  images{ii} = rgb2gray(currentimage);
  #figure, imshow(images{ii})
end


sigma = sqrt(2)/2;
gauss1 = fspecial('gaussian', round([10*sigma 10*sigma]), sigma);
sigma = sqrt(2);
gauss2 = fspecial('gaussian', round([10*sigma 10*sigma]), sigma);
blur1 = imfilter(images{1}, gauss1, 'replicate', 'same');
blur2 = imfilter(images{1}, gauss2, 'replicate', 'same');
dog2 = blur1 - blur2;

subplot(1,2,1),
imshow(images{1})

subplot(1,2,2),
imshow(dog2)

I = images{1};

%% Filter using DoG
stepsPerOctave = 5;
octaves = 2;
#mult = nthroot(2,stepsPerOctave);
mult = sqrt(2);

% Create blurry images
sigma = sqrt(2)/2;
kernelSize = [ceil(6*sqrt(2)*sigma),ceil(6*sqrt(2)*sigma)]

tx = ty = linspace (-2, 2, 5)';
[xx, yy] = meshgrid (tx, ty);

J = images{1};

for oct=0:octaves-1
  
  I = J;

  for k = 1:stepsPerOctave
      
      step = oct*(stepsPerOctave)+k;
      disp(['sigma_cal is ' num2str(sigma)]), disp(['iter is ' num2str(step)]);
      
      gauss = fspecial('gaussian', kernelSize, sigma);
      blur{step} = imfilter(I, gauss, 'replicate', 'same');
      sigma = sigma * mult;
      
      #{
      figure, subplot(1,2,1), imagesc(blur{oct*stepsPerOctave+k}); colorbar; 
      title(['Gaussian ' num2str(oct*stepsPerOctave+k) num2str(sigma)]); #pause;
      
      subplot(1,2,2)
      mesh (tx, ty, gauss);
      xlabel ("tx");
      ylabel ("ty");
      zlabel ("value");
      title ("3-D gauss plot");
      #}
  end
  sigma = sigma / power(mult,3);
  gauss = fspecial('gaussian', kernelSize, sigma);
  disp(['sigma_cal is ' num2str(sigma)]),
  J = imfilter(I, gauss, 'replicate', 'same');
  J =   imresize(J, 0.5);
  
end

sigma = sqrt(2)/2;
% Create DoG
for oct = 0:octaves-1
  for k = 1:stepsPerOctave-1
      
      step = oct*(stepsPerOctave)+k
      disp(['iter is ' num2str(step)]);
      dog{step} = blur{step+1}-blur{step};
      
      figure, imagesc(dog{step}); colorbar; 
      title(['DoG ' num2str(step) '   ' num2str(sigma)]);
      
      disp(['iter ' num2str(step) ' over']);
      
      sigma = sigma * mult;
  end
  sigma = sigma / power(mult,2);
end


