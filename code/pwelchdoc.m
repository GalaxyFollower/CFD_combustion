%Obtain Welch's overlapped segment averaging PSD estimate of the preceding signal. 
%Use a segment length of 500 samples with 300 overlapped samples. 
%Use 500 DFT points so that 100 Hz falls directly on a DFT bin. DISCRETE
%FOURIER TRANSFORM
%Input the sample rate to output a vector of frequencies in Hz. Plot the result.

fs = 1000;
t = 0:1/fs:5-1/fs;
x = cos(2*pi*100*t);

[pxx,f] = pwelch(x,500,300,500,fs);

plot(f,10*log10(pxx))
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

[pks,locs] = findpeaks(pxx);
f(locs)%has to be 100 here (the frequency)