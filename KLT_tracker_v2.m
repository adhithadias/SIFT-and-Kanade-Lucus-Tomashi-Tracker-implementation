#function KLT_tracker_v2

close all

folder = '.\images4';
im_array = 0:7;
starting_frame = 845;
num_of_images = numel(im_array);
images = cell(num_of_images,1);
for k = im_array,
    #images{k+1} = imread(fullfile(folder, ['hotel.seq' num2str(k) '.png']));
    images{k+1} = imread(fullfile(folder, ['frame' num2str(starting_frame+k) '.jpg']));
    images{k+1} = rgb2gray(images{k+1});
    images{k+1} = im2single(images{k+1});
end

%%
% detect keypoints
% calculate harris corners
tau = 0.06;
num_of_points = 500;
i1 = images{1};

gauss_filter = fspecial('gaussian', [7 7], 1);
imblur = imfilter(i1, gauss_filter);
[Ix, Iy] = gradient(imblur);
Ixx = imfilter(Ix.*Ix, gauss_filter);   % second derivative w.r.t. x
Iyy = imfilter(Iy.*Iy, gauss_filter);   % second derivative w.r.t. y
Ixy = imfilter(Ix.*Iy, gauss_filter);   % second derivative w.r.t. x and y
har = Ixx.*Iyy - Ixy.*Ixy - tau*(Ixx+Iyy).^2;   % Harris criterion

% check that center is strict max
maxv = ordfilt2(har, 49, ones(7)); % maximum filter
maxv2 = ordfilt2(har, 48, ones(7));   % second max after ordering in ascending order
ind = find(maxv==har & maxv~=maxv2);  % check if it is a local maximum

% get top N points
[sv, sind] = sort(har(ind), 'descend');
sind = ind(sind);
[pty, ptx] = ind2sub(size(i1), sind(1:min(num_of_points, numel(sind))));
%[pty, ptx] = find(maxv==har & maxv>0);
figure(1), imagesc(i1), hold off, axis image, colormap gray
hold on, plot(ptx, pty, 'r.'); drawnow;

% Track points
ws = 7;
% initialize
track_x = zeros(num_of_points, num_of_images);
track_y = zeros(num_of_points, num_of_images);
track_x(:, 1) = ptx(:);  % initially selected points
track_y(:, 1) = pty(:);  % initially selected points


for t = 1:num_of_images-1
    [track_x(:, t+1), track_y(:, t+1)] = ...
            getNextPoints(track_x(:, t), track_y(:, t), images{t}, images{t+1}, ws);
end



#{
#########################################################################
#########     KLT Tracker algorithm Implementation
#########################################################################

im1 = images{1};
im2 = images{2};
x = track_x(:,1);
y = track_y(:,1);

x2 = x;
y2 = y;

hw = floor(ws/2);
[winx, winy] = meshgrid(-hw:hw, -hw:hw);
winx = winx(:);
winy = winy(:);

gaus_filt = fspecial('gaussian', [7 7], 1);
[Ix, Iy] = gradient(imfilter(im1, gaus_filt));

% do interpolation in vectorized form
xw = repmat(x, [1 numel(winx)])'+repmat(winx', [numel(x) 1])';
disp( ['xw ', num2str(size(xw))] );
yw = repmat(y, [1 numel(winy)])'+repmat(winy', [numel(y) 1])';

valid = ~isnan(x) & x>=hw & y>=hw & x+hw<=size(im1,2) & y+hw<=size(im1,1);
disp(class(valid))
disp( ['valid size ' num2str(size(valid))] );
patch1_all = zeros(size(xw));
Ix_all = zeros(size(xw));
Iy_all = zeros(size(xw));
patch1_all(:, valid) = interp2(im1, xw(:, valid), yw(:, valid), 'cubic');
Ix_all(:, valid) = interp2(Ix, xw(:, valid), yw(:, valid), 'cubic');
Iy_all(:, valid) = interp2(Iy, xw(:, valid), yw(:, valid), 'cubic');

for iter = 1:5

    % need to make sure point hasn't drifted too close to border
    valid2 = valid & x2>=hw & y2>=hw & x2+hw<=size(im2,2) & y2+hw<=size(im2,1);
    x2(~valid2) = nan;  # replace valid2 zeros with not a number
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
        x2(p)=x2(p)+d(1)
        y2(p)=y2(p)+d(2);        
       
    end
end

track_x(:,2) = x2;
track_y(:,2) = y2;

###############################################################################

#}

figure(2), imagesc(images{1}), hold off, axis image, colormap gray
hold on, plot(track_x', track_y', 'r')
figure(3), imagesc(images{end}), hold off, axis image, colormap gray
hold on, plot(track_x', track_y', 'r')