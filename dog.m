clear all
close all

sigma = 3;
t = [-10:0.05:10];
k = 0.4;

gauss = gaussian(t,sigma);
first_derivative_of_gauss = -t/power(sigma,2).*gauss;

hold on
plot(t,gauss);
plot(t,first_derivative_of_gauss);
hold off

second_derivative_of_gauss = (t.*t-power(sigma,2)).*gauss/power(sigma,4);

figure,
plot(t,second_derivative_of_gauss);
title("Difference of Gaussian for different k values");
xlabel("t");
ylabel("value");
grid on

#hold on

for i=1:1:10
  k = 1.2 + i*0.2;
  if(k!=1)
    disp ("The value of iteration is:"), disp (k),
    diff_of_gauss = difference_of_gaussian(t,sigma,k);
    subplot(2,5,i);
    title(['Difference of Gaussian for k = ' num2str(k)]);
    xlabel("t");
    ylabel("value");
    grid on
    hold on
    plot(t,diff_of_gauss); 
    plot(t,second_derivative_of_gauss);
    hold off
  else
    disp ("The value of iteration is:"), disp (k),
    subplot(2,5,i);
    title(['Difference of Gaussian for k = ' num2str(k)]);
    xlabel("t");
    ylabel("value");
    grid on
    hold on
    plot(t,diff_of_gauss); 
    #plot(t,second_derivative_of_gauss);
    hold off
  endif
end

#{
disp ("The value of iteration is:"), disp (k),
diff_of_gauss = difference_of_gaussian(t,sigma,k);
plot(t,diff_of_gauss);
#}
  
#hold off

