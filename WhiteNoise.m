function [ U ] = WhiteNoise(N,fraction_f_nyq,mu,variance)
%uit N samples kan je maximaal sinussen met nyquistfreq=N/2 voorstellen
%de band laat frequenties door tot fraction_f_nyq*f_nyq
Type='RGS';
Band=[0 fraction_f_nyq];
ulow=mu-sqrt(variance);
uhigh=mu+sqrt(variance);
Levels=[ulow uhigh];
U=idinput(N,Type,Band,Levels);
% figure(12);
% plot(U);
%  E=pwelch(U);
%  figure(2);
%  plot(E);
%  F=fft(U);
%  figure(3);
%  semilogy(abs(F));
end

