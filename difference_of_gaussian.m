function dog = difference_of_gaussian(t,sigma,k);
  
  dog = (gaussian(t,k*sigma)-gaussian(t,sigma)) / ((k-1)*power(sigma,1));
  
endfunction