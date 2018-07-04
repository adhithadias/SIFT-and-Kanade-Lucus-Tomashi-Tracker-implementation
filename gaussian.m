function gauss = gaussian(t,sigma)
  
  gauss = 1/sqrt(2*pi)/sigma * exp(-(t.*t)/(2*power(sigma,2)));

endfunction
