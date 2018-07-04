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


sigma = 0.5;
gauss1 = fspecial('gaussian', round([10*sigma 10*sigma]), sigma);
sigma = 1;
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
mult = nthroot(2,stepsPerOctave);

% Create blurry images
sigma = 0.5;
kernelSize = [ceil(6*sqrt(2)*sigma),ceil(6*sqrt(2)*sigma)]

tx = ty = linspace (-2, 2, 5)';
[xx, yy] = meshgrid (tx, ty);

#{
for k = 1:octaves*stepsPerOctave+1
    disp(['Sigma is ' num2str(sigma)]);
    gauss = fspecial('gaussian', kernelSize, sigma);
    blur(:,:,k) = imfilter(I, gauss, 'replicate', 'same');
    #figure, imagesc(blur(:,:,k)); colorbar; title(['Gaussian ' num2str(k)]); #pause;
    figure, subplot(1,2,1), imshow(blur(:,:,k)), subplot(2,1,2),
    
    mesh (tx, ty, gauss);
    xlabel ("tx");
    ylabel ("ty");
    zlabel ("value");
    title ("3-D gauss plot");
    
    sigma = sigma * mult;
end
#}

J = images{1};

for oct=0:octaves-1
  
  I = J;

  for k = 1:stepsPerOctave+1
      
      step = oct*(stepsPerOctave+1)+k;
      disp(['sigma_cal is ' num2str(sigma)]), disp(['iter is ' num2str(step)]);
      
      gauss = fspecial('gaussian', kernelSize, sigma);
      blur{step} = imfilter(I, gauss, 'replicate', 'same');
      #blur(:,:,k,oct) = imfilter(I, gauss, 'replicate', 'same');
      sigma = sigma * mult;
      
      #{
      figure, subplot(1,2,1), imagesc(blur{oct*stepsPerOctave+k}); colorbar; title(['Gaussian ' num2str(oct*stepsPerOctave+k)]); #pause;
      #figure, subplot(1,2,1), imagesc(blur(:,:,k)); colorbar; title(['Gaussian ' num2str(k)]); #pause;
      #figure, subplot(1,2,1), imshow(blur(:,:,k)),
      
      subplot(1,2,2)
      mesh (tx, ty, gauss);
      xlabel ("tx");
      ylabel ("ty");
      zlabel ("value");
      title ("3-D gauss plot");
      #}
  end
  disp(['sigma_cal is ' num2str(sigma)]),
  #gauss = fspecial('gaussian', kernelSize, sigma);
  J = imfilter(I, gauss, 'replicate', 'same');
  J =   imresize(J, 0.5);
  
end

#{
% Create DoG
for k = 1:stepsPerOctave
    #dog(:,:,k,1) = blur(:,:,k+1) - blur(:,:,k);
    #imagesc(dog(:,:,k)); colorbar; title(['DoG ' num2str(k)]); pause;
    dog{k} = blur{k+1}-blur{k};
    imagesc(dog(:,:,k)); colorbar; title(['DoG ' num2str(k)]); pause;
end
#}

% Create DoG
for oct = 0:octaves-1
  for k = 1:stepsPerOctave
      #dog(:,:,k,1) = blur(:,:,k+1) - blur(:,:,k);
      #imagesc(dog(:,:,k)); colorbar; title(['DoG ' num2str(k)]); pause;
      
      step = oct*(stepsPerOctave+1)+k
      disp(['iter is ' num2str(step)]);
      dog{step} = blur{step+1}-blur{step};
      figure, imagesc(dog{step}); colorbar; title(['DoG ' num2str(step)]); #pause;
      disp(['iter ' num2str(step) ' over']);
  end
end




