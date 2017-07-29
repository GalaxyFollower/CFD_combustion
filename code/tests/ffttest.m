close all
clear all

freqspec = zeros(1,10000);
freqspec(101)=1000;

x = ifft(freqspec);
plot(real(x));
figure;

timelength = 0.1;
timescale = linspace(0,timelength,length(freqspec));
delta_t = timescale(10)-timescale(9);
y = fft(x);
plot(1/delta_t*timescale',real(y'));