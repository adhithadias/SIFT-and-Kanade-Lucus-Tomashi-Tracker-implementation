function KLT_tracker

folder = '.\images';
im = readImages(folder, 0:20);

tau = 0.06;
N = 250;
[pt_x, pt_y] = detectKeypoints(im{1}, tau, N);

ws = 7;
[track_x, track_y] = trackPoints(pt_x, pt_y, im, ws);
  
figure(2), imagesc(im{1}), hold off, axis image, colormap gray
hold on, plot(track_x', track_y', 'r')
figure(3), imagesc(im{end}), hold off, axis image, colormap gray
hold on, plot(track_x', track_y', 'r')




function im = readImages(folder, nums)
  im = cell(numel(nums),1);
  t = 0;
  for k = nums,
      t = t+1;
      im{t} = imread(fullfile(folder, ['hotel.seq' num2str(k) '.png']));
      im{t} = im2single(im{t});
  end


function [ptx, pty] = detectKeypoints(im, tau, N)

  % get harris function
  gfil = fspecial('gaussian', [7 7], 1);
  imblur = imfilter(im, gfil);
  [Ix, Iy] = gradient(imblur);
  Ixx = imfilter(Ix.*Ix, gfil);
  Iyy = imfilter(Iy.*Iy, gfil);
  Ixy = imfilter(Ix.*Iy, gfil);
  har = Ixx.*Iyy - Ixy.*Ixy - tau*(Ixx+Iyy).^2;

  % check that center is strict max
  maxv = ordfilt2(har, 49, ones(7)); % sorts values in each window
  maxv2 = ordfilt2(har, 48, ones(7));
  ind = find(maxv==har & maxv~=maxv2); 

  % get top N points
  [sv, sind] = sort(har(ind), 'descend');
  sind = ind(sind);
  [pty, ptx] = ind2sub(size(im), sind(1:min(N, numel(sind))));
  %[pty, ptx] = find(maxv==har & maxv>0);
  figure(1), imagesc(im), hold off, axis image, colormap gray
  hold on, plot(ptx, pty, 'r.'); drawnow;


function [track_x, track_y] = trackPoints(pt_x, pt_y, im, ws)

  % initialize
  N = numel(pt_x);
  nim = numel(im);
  track_x = zeros(N, nim);
  track_y = zeros(N, nim);
  track_x(:, 1) = pt_x(:);
  track_y(:, 1) = pt_y(:);

  for t = 1:nim-1
      [track_x(:, t+1), track_y(:, t+1)] = ...
              getNextPoints(track_x(:, t), track_y(:, t), im{t}, im{t+1}, ws);
  end


function [x2, y2] = getNextPoints(x, y, im1, im2, ws)

  x2 = x;
  y2 = y;

  hw = floor(ws/2);
  [winx, winy] = meshgrid(-hw:hw, -hw:hw);
  winx = winx(:);
  winy = winy(:);

  gausfil = fspecial('gaussian', [7 7], 1);
  [gx, gy] = gradient(imfilter(im1, gausfil));

  % do interpolation in vectorized form
  xw = repmat(x, [1 numel(winx)])'+repmat(winx', [numel(x) 1])';
  yw = repmat(y, [1 numel(winy)])'+repmat(winy', [numel(y) 1])';

  valid = ~isnan(x) & x>=hw & y>=hw & x+hw<=size(im1,2) & y+hw<=size(im1,1);
  patch1_all = zeros(size(xw));
  Ix_all = zeros(size(xw));
  Iy_all = zeros(size(xw));
  patch1_all(:, valid) = interp2(im1, xw(:, valid), yw(:, valid), 'cubic');
  Ix_all(:, valid) = interp2(gx, xw(:, valid), yw(:, valid), 'cubic');
  Iy_all(:, valid) = interp2(gy, xw(:, valid), yw(:, valid), 'cubic');

  for iter = 1:5

      % need to make sure point hasn't drifted too close to border
      valid2 = valid & x2>=hw & y2>=hw & x2+hw<=size(im2,2) & y2+hw<=size(im2,1);
      x2(~valid2) = nan;
      y2(~valid2) = nan;

      x2w = repmat(x2, [1 numel(winx)])'+repmat(winx', [numel(x) 1])';
      y2w = repmat(y2, [1 numel(winy)])'+repmat(winy', [numel(y) 1])';  
      patch2_all = zeros(size(patch1_all));
      patch2_all(:, valid2) = interp2(im2, x2w(:, valid2), y2w(:, valid2), 'cubic');    
      
      % for each point, step towards the matching location
      for p = 1:numel(x)

          if ~valid2(p)
              continue;
          end
          
          patch1 = patch1_all(:, p);
          Ix = Ix_all(:, p);
          Iy = Iy_all(:, p);
          patch2 = patch2_all(:, p);       

          It = patch2-patch1;
          A = [sum(Ix.*Ix) sum(Ix.*Iy) ;
               sum(Ix.*Iy) sum(Iy.*Iy)];
          b = -[sum(Ix.*It) ; sum(Iy.*It)];
          d = A\b; 
          x2(p)=x2(p)+d(1);
          y2(p)=y2(p)+d(2);        

      end
  end