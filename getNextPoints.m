function [x2, y2] = getNextPoints(x, y, im1, im2, ws)

  x2 = x;
  y2 = y;

  hw = floor(ws/2);
  [winx, winy] = meshgrid(-hw:hw, -hw:hw);
  winx = winx(:);
  winy = winy(:);

  gaus_filt = fspecial('gaussian', [7 7], 1);
  [gx, gy] = gradient(imfilter(im1, gaus_filt));

  % do interpolation in vectorized form
  xw = repmat(x, [1 numel(winx)])'+repmat(winx', [numel(x) 1])';
  disp( ['xw ', num2str(size(xw))] );
  yw = repmat(y, [1 numel(winy)])'+repmat(winy', [numel(y) 1])';

  valid = ~isnan(x) & x>=hw & y>=hw & x+hw<=size(im1,2) & y+hw<=size(im1,1);
  disp(class(valid))
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
  
  