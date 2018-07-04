close all

imagefiles = dir('graf/*.ppm');      
nfiles = length(imagefiles);    % Number of files found

for ii=1:nfiles
  disp( ['Reading image ' num2str(ii)] );
  currentimage = imread(['graf/' imagefiles(ii).name]);
  images{ii} = rgb2gray(currentimage);
  images{ii} = im2single(images{ii});
  #images{ii} = currentimage;
  #figure, imshow(images{ii})
end

imshow(images{1})

I = images{1};
figure(1), imagesc(I), hold off, axis image, colormap gray

%% Filter using DoG
stepsPerOctave = 5;
octaves = 2;
# mult = nthroot(2,stepsPerOctave);
mult = sqrt(2);

% Create blurry images
sigma = sqrt(2)/2;
kernelSize = [ceil(6*sqrt(2)*sigma),ceil(6*sqrt(2)*sigma)]

tx = ty = linspace (-2, 2, 5)';
[xx, yy] = meshgrid (tx, ty);

J = images{1};

for oct=0:octaves-1
  
  I = J;
  disp(['Sigma at the start of octave ' num2str(oct+1) ' is ' num2str(sigma)]);

  for k = 1:stepsPerOctave
      
      step = oct*(stepsPerOctave)+k;
      disp(['With sigma = ' num2str(sigma) ' creating blurred image step ' num2str(step) ' octave: ' num2str(oct+1) ' and level: ' num2str(k)]);
      
      gauss = fspecial('gaussian', kernelSize, sigma);
      blur{step} = imfilter(I, gauss, 'replicate', 'same');
      sigma = sigma * mult;
      
      tau = 0.06;
      num_of_points = 500;
      [Ix, Iy] = gradient(blur{step});
      Ixx = imfilter(Ix.*Ix, gauss);   % second derivative w.r.t. x
      Iyy = imfilter(Iy.*Iy, gauss);   % second derivative w.r.t. y
      Ixy = imfilter(Ix.*Iy, gauss);   % second derivative w.r.t. x and y
      har = Ixx.*Iyy - Ixy.*Ixy - tau*(Ixx+Iyy).^2;
      harris{step} = har;
      
      % check that center is strict max
      maxv = ordfilt2(har, 49, ones(7)); % maximum filter
      maxv2 = ordfilt2(har, 48, ones(7));   % second max after ordering in ascending order
      ind = find(maxv==har & maxv~=maxv2);  % check if it is a local maximum
      
      indices{step} = ind;
      
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
  J = imfilter(I, gauss, 'replicate', 'same');
  J = imresize(J, 0.5);
  
end

sigma = sqrt(2)/2;
% Create DoG
for oct = 0:octaves-1
  for k = 1:stepsPerOctave-1
      
      step = oct*(stepsPerOctave)+k;
      disp(['DOG calculation step: ' num2str(step) ' Octave: ' num2str(oct+1) ' Level: ' num2str(k)]);
      dog{step} = blur{step+1}-blur{step};
      
      [m, n] = size(dog{step});
      ind = indices{step};
      [pty, ptx] = ind2sub(size(dog{step}), ind);
      figure, imagesc(dog{step}), hold off, axis image, colormap gray
      hold on, plot(ptx, pty, 'r.'); drawnow;
      #figure, imshow(dog{step}); 
      
      title(['DoG calculation step ' num2str(step) ' with sigma ' num2str(sigma)]);
      
      sigma = sigma * mult;
  end
  sigma = sigma / power(mult,2);
end

hold off
% Detect keypoints using the scale-space extrema detection followed by Harris Corner Criterion
sigma = sqrt(2)/2;
D = [];   
D_new = [];
D_all = [];
D_harris = [];                            # keypoints
#for k = 1:stepsPerOctave-1

sigma = sqrt(2)/2;

for k = 1:stepsPerOctave-1
  D = [];   
  D_new = [];
  disp(['Keypoint calculation iteration number: ' num2str(k)]);
  DOG = dog{k};
  [m,n] = size(DOG);
  ind = indices{k};

  if(k == 1)
    lower_level = dog{k+1};
    middle_level = dog{k};
    upper_level = dog{k+1};          
  elseif(k == stepsPerOctave-1)
    lower_level = dog{k-1};
    middle_level = dog{k};
    upper_level = dog{k-1};
  else
    lower_level = dog{k-1};
    middle_level = dog{k};
    upper_level = dog{k+1};    
  endif

  for i=2:m-1
    for j=2:n-1
  
      mat(:,:,1) = lower_level(i-1:i+1,j-1:j+1);
      mat(:,:,2) = middle_level(i-1:i+1,j-1:j+1);               
      mat(:,:,3) = upper_level(i-1:i+1,j-1:j+1); 
  
      val = mat(2,2,2);             # replace middle value with a neighboring value
      mat(2,2,2) = mat(1,1,1);
      
      mat_min = min(min(min(mat)));
      mat_max = max(max(max(mat)));
      
      if (val>mat_max || val<mat_min)
          d = [i; j; sigma];
          D = [D d];
          D_all = [D_all d];
         
          if ( ismember((i-1)*n + j, ind) && i-3>0 && i+3<m && j-3>0 && j+3<n)
            disp([num2str(i) '  ' num2str(j) '  ' num2str(i*m + j)]);
            D_new = [D_new d];  
            D_harris = [D_harris d];
          endif 
      end
      
    endfor  
  endfor

  disp( [ 'Keypoints by the end of iteration: ' num2str(size(D_harris)(2)) '  sigma: ' num2str(sigma) ] );
  sigma = sigma * mult;
  
  keypD{k} = D;
  keypD_harris{k} = D_new;
  
end

disp( ['Size of keypoints is ' num2str(size(D_harris)(2))] );

# Outlier rejection using Harris criterion
sigma = sqrt(2)/2;

disp( ['Size of new keypoints is ' num2str(size(D_new)(2))] );

figure, imagesc(images{1}), hold off, axis image, colormap gray
hold on, scatter(D_harris(2,:)',D_harris(1,:)',10*D_harris(3,:)','r'); drawnow;

D_new = [];

for p=1:size(D_harris)(2)

    hist = zeros(1, 36);        # 36 positions in the histogram
    i = D_harris(1,p);   j = D_harris(2,p);   sig = D_harris(3,p);
    #D_harris(:,p)
    sigma = sqrt(2)/2;
    
    I = images{1};
    for k = 1:stepsPerOctave-1
      if(sig == sigma)
        I = blur{k};
        break;
      endif
      sigma = sigma * mult;
    endfor
    
    hw = 2;
    gauss = fspecial('gaussian', [2*hw+1,2*hw+1], sig*1.5);

    m = 0;  theta = 0;
    for k=-hw:hw
        for l=-hw:hw
            #[m, theta] = point_grad(I, i+k, j+l);
            x = i+k;    y = j+l;
            
            m = sqrt( (I(x,y-1)-I(x,y+1))^2 + (I(x-1, y)-I(x+1, y))^2 );
            theta = atan( (I(x+1, y)-I(x-1, y))/(I(x, y+1)-I(x, y-1)) );
            if theta<0
                theta = 2*pi+theta;
            end
            theta = theta*180/pi;
            
            idx = idivide(theta, 10);
            hist(idx+1) = hist(idx+1) + m * gauss(k+hw+1, l+hw+1);
        end
    end
    %     hist
    [max_grad,orientation] = max(hist);
    
    D_harris(4,p) = (orientation-1)*10;
    D_harris(5,p) = max_grad*1000*sig+sig;
    
    hist(orientation) = -Inf;
    [second_max,orient] = max(hist);
    
    if(second_max > max_grad*0.8)
      d = [i; j; sig; (orient-1)*10; second_max*1000*sig+sig];
      disp([" A SECOND ORIENTATAION FOUND " num2str((orient-1)*10)]);
      D_new = [D_new d];
    endif

    
    disp([num2str(D_harris(1,p)) '  ' num2str(D_harris(2,p)) '  ' num2str(D_harris(3,p)) '  ' num2str(D_harris(4,p)) '  ' num2str(D_harris(5,p))]);
endfor

D_harris = [D_harris D_new];
disp( [num2str(size(D_new)(2)) ' number of SECOND ORIENTATAIONS ADDED'] );

x1 = []; x2 = []; y1 = []; y2 = [];
arrow = [];
for p=1:size(D_harris)(2)
  i = D_harris(1,p);   j = D_harris(2,p);   sig = D_harris(3,p); orient = D_harris(4,p)/180*pi;
  mag = D_harris(5,p);
  
  i2 = i - mag / 5 * cos(orient); j2 = j + mag / 5 * sin(orient);
  ar = [j, i, j2, i2];
  arrow = [arrow; ar];
endfor

figure, imagesc(images{1}), hold off, axis image, colormap gray
hold on, scatter(D_harris(2,:)',D_harris(1,:)',D_harris(5,:)','y');
drawArrow(arrow, 2, 0.7, 5, 0.5); drawnow;