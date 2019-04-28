clear;
close all;

Nfft = 2^13;
fs = 10e3;
f0 = fs/600;
s = cos(2 * pi * f0 * (0:1/fs:1-1/fs));
s = int32(s*(2^10-1));
s = complex(s);
WriteDataIQ('~/projectsDSP/commonpocket/Debug/CSC_cmake/Bin/cosf0.bin', s, 'int32');
spec = abs(fft(s, Nfft));

% figure(1); plot(real(s)); grid on;
figure(2); plot(spec); grid on;