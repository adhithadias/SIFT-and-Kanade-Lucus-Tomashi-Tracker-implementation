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
            #theta = theta*180/pi;
            
            idx = idivide(theta, 10);
            hist(idx+1) = hist(idx+1) + m * gauss(k+hw+1, l+hw+1);
        end
    end
    %     hist
    [max_grad,orientation] = max(hist);
    
    D_harris(4,p) = (orientation-1)*10;
    D_harris(5,p) = max_grad*255*D_harris(2,p)+D_harris(2,p);
    
    hist(orientation) = -Inf;
    [second_max,orient] = max(hist);
    
    if(second_max > max_grad*0.8)
      d = [i; j; sig; (orient-1)*10; second_max*255*sig+sig];
      disp([" A SECOND ORIENTATAION FOUND " num2str((orient-1)*10)]);
      D_new = [D_new d];
    endif

    
    disp([num2str(D_harris(1,p)) '  ' num2str(D_harris(2,p)) '  ' num2str(D_harris(3,p)) '  ' num2str(D_harris(4,p)) '  ' num2str(D_harris(5,p))]);
endfor

D_harris = [D_harris D_new];
disp( [num2str(size(D_new)(2)) ' number of SECOND ORIENTATAIONS ADDED'] );

x1 = []; x2 = []; y1 = []; y2 = [];
arrow = []
for p=1:size(D_harris)(2)
  i = D_harris(1,p);   j = D_harris(2,p);   sig = D_harris(3,p); orient = D_harris(4,p);
  mag = D_harris(4,p);
  
  i2 = i - mag * cos(orient); j2 = j + mag * sin(orient);
  arrow = [arrow; j, i, j2, i2];
endfor

figure, imagesc(images{1}), hold off, axis image, colormap gray
#hold on, scatter(D_harris(2,:)',D_harris(1,:)',D_harris(5,:)','r'); drawnow;
hold on, drawArrow(arrow, 1, 0.1); drawnow;