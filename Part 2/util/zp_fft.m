function [Z_freq] = zp_fft(Z_time,k)
% Make sure K is greater than half the length of Z_time
len = length(Z_time);
x1 = Z_time(1:k);
x2 = Z_time(k+1:len);
x2 = [x2 zeros(1,k-length(x2))];
x = x1 + x2;
Z_freq = fft(x);
end
