addpath(genpath('../sc_common'));
addpath(genpath('../matlab_utils'));
addpath(genpath('../sc_amra'));
clear;
close all;

fsHz = 80000;

%%
% sDmr = ReadDataIQ('signals/sDmr19200.bin', 'int32');
% fsDmrHz = 19200;
% [p, q] = rat(max(fsHz, fsDmrHz) / min(fsHz, fsDmrHz));
% sDmrRsmpl = resample(sDmr, p, q);
% WriteDataIQ(strcat('signals/sDmr', num2str(fsHz), '.bin'), sDmrRsmpl, 'int32');

% %%
% T = 0.1;
% f0 = 400;
% snr = 20;
% ampl = (2^10-1);
% sCos = cos(2 * pi * f0 .* (0 : 1/fsHz : T-1/fsHz));
% 
% % sDSBf0_80000 = complex(sCos);
% % sDSBf0_80000_snr = awgn(sDSBf0_80000, snr);
% % WriteDataIQ(strcat('signals/sDSB', num2str(f0), '_', num2str(fsHz), '_snr', num2str(snr), 'dB'), ...
% %     ampl * sDSBf0_80000_snr, 'int32');
% % figure(1);
% % plot(real(sDSBf0_80000_snr));
% % grid on;
% % 
% sAMf0_80000 = complex((1 + 0.5 * sCos));
% sAMf0_80000_snr = awgn(sAMf0_80000, snr);
% % WriteDataIQ(strcat('signals/sAM', num2str(f0), '_', num2str(fsHz), '_snr', num2str(snr), 'dB'), ... 
% %     ampl * sAMf0_80000_snr, 'int32');
% % figure(2);
% % plot(real(sAMf0_80000_snr));
% % grid on;
% % 
% % fDevHz = fsHz / 16;
% % fcHz = fDevHz * 2;
% % xFM = fmmod(sCos, fcHz, fsHz, fDevHz);
% % sFMf0_80000 = hilbert(xFM) .* exp(-1i * 2*pi*fcHz * (0:length(xFM)-1)/fsHz);
% % sFMf0_80000_snr = awgn(sFMf0_80000, snr);
% % WriteDataIQ(strcat('signals/sFM', num2str(fDevHz), '_', num2str(fsHz), '_snr', num2str(snr), 'dB'), ... 
% %     ampl * sFMf0_80000_snr, 'int32');
% % 
% % sZeros = complex(zeros(1, length(sCos)));
% % sZeros_snr = awgn(sZeros, snr);
% % WriteDataIQ(strcat('signals/sZeros', '_', num2str(fsHz), '_snr', num2str(snr), 'dB'), ... 
% %     ampl * sZeros_snr, 'int32');
% 
% % figure(3);
% % plot(real(sFMf0_80000_snr));
% % hold on;
% % % plot(imag(sFMf0_80000_snr));
% % % hold on;
% % plot(sCos);
% % grid on;
% 
% % Nfft = 2^nextpow2(length(xFM));
% % % Nfft = length(xFM);
% % freqs = (-Nfft/2 : (Nfft-1)/2) / Nfft * fsHz;
% % figure(4);
% % plot(freqs, mag2db(fftshift(abs(fft(sFMf0_80000_snr, Nfft)))));
% % % plot(freqs, (fftshift(abs(fft(xFM, Nfft)))));
% % grid on;

%%
% env = ReadDataIQ('signals/pocket/tmp/defaultSource_413749999_413753126_1ch_40id.bin', 'int32');
% env = ReadDataIQ('signals/pocket/tmp/defaultSource_401249999_401253126_1ch_44id.bin', 'int32');
% env = ReadDataIQ('signals/pocket/raw/Raw__165_90100000_200000__468750__18_45_23_961', 'int16');
% env = ReadDataIQ('signals/pocket/raw/Raw__165_66300000_40000__78125__19_15_48_697', 'int16');
% envelope = sAMf0_80000_snr;
fsHz = 13333.333;
% fsHz = 15625;

folder = 'signals/pocket/';
%--am--
% file = 'am/Raw__41_153001562_10000__13333__03_02_29_337';
% file = 'am/Raw__41_153001562_10000__13333__03_03_03_560';
% file = 'am/Raw__40_157000000_10000__13333_';
% file = 'am/Raw__40_157001562_10000__13333_';
%--fm--
file = 'fm/Raw__41_436001562_10000__13333__03_00_41_826';
% file = 'fm/Raw__165_435001000_100000__234375__20_45_38_949';

% env = ReadDataIQ('~/VirtusData/Default/VirtusRecordedRawData/2019_02_13/Raw__165_157000000_10000__15625__03_00_00_000', 'int16')';

env = ReadDataIQ([folder, file], 'int16')';
% env = env(1:end/4);
% pos = 1;
% env = env(pos : pos + 4096);
% env = env(110540:134000);
t = (0 : length(env)-1) / fsHz;

%%
% lenF = 2^11;
% s = zeros(1, ceil(length(env) / lenF) * lenF);
% s(1:length(env)) = env;
% s = reshape(s, lenF, length(s)/lenF)';
% spec = fftshift(abs(fft(s, lenF, 2)), 2);
% figure(7);
% mesh(mag2db(spec));
% grid on;

%%
% env = env(1500:end);
% WriteDataIQ([folder, file, '-iq32'], env, 'int32');

%%
m = mean(abs(env));
iNW = 1;
for i = 2 : length(env)
    if (abs(env(i)) > m/4 && abs(env(i-1)) > m/4) 
        phDiff(iNW) = angle(env(i) * conj(env(i-1)));
        iNW = iNW + 1;
    end
end
% phDiff = angle(env(2:end) .* conj(env(1:end-1)));
phDiffMean = mean(phDiff);
phDiffStd = std(phDiff);
phDiffAbsStd = std(abs(phDiff));
fOffset = phDiffMean * fsHz / (2*pi)
% env = env .* exp(1i * (-2*pi*fOffset) * t);

freqStdHz = phDiffStd * fsHz / (2*pi)


%%
% s = amdemod(real(env), abs(fOffset), fsHz); 

%%
fs = 13333;
N = 2^nextpow2(length(env));
freqs = (-N/2 : (N-1)/2) * fs / N;
figure(1);
subplot(2,1,1); plot(t, real(env)); hold on; plot(t, imag(env)); grid on; xlabel('time, s');
subplot(2,1,2); plot(freqs, mag2db(fftshift(abs(fft(env, N))))); grid on; xlabel('frequency, Hz');
set(gcf, 'color', 'w'); set(groot, 'DefaultAxesFontSize', 18);

figure(2);
plot(freqs, mag2db(fftshift(abs(fft(env/max(abs(env)) / (N/2), N))))); 
grid on; 
xlabel('Frequency, Hz');
ylabel('|S(f)|, dB')

% figure(8);
% subplot(2,1,1);
% plot(real(env)); hold on; plot(imag(env)); grid on;
% subplot(2,1,2);
% plot(phDiff);
% grid on;

return;

%%
thresholds.ampl = 1;
thresholds.P = 0.2;
thresholds.gammaMax = 0.6e-3;
thresholds.sigmaAP = pi/4;
thresholds.sigmaDP = pi/2 * 0.9;

kf = KeyFeatures(env, thresholds.ampl)
decision = AMRA1(kf, thresholds)

% figure(2);
% plot(abs(env));
% grid on

env1 = env(abs(env) > 3);
phDiff = env(2:end) .* conj(env(1:end-1));
phDiff = angle(phDiff(abs(env(2:end)) > 50));
figure(5);
plot(t(1:length(phDiff)), phDiff * fsHz / (2*pi));
grid on;
set(gcf, 'color', 'w'); set(groot, 'DefaultAxesFontSize', 18);
xlabel('Time, s');
ylabel('Inst. frequency, Hz')

phase = angle(env);
phase = unwrap(phase);
phase = phase - (2*pi*fOffset*t);
phStd = std(phase)
figure(6);
plot(phase);
hold on;
plot(phase - mean(phase));
grid on;

mult = env(2:end) .* conj(env(1:end-1));
f1 = abs(sum(mult));
f2 = sum(abs(mult));
m = f1 / f2

% fs1 = 468750;
% fs2 = 48000;
% factor = fs2 / fs1;
% [p, q] = rat(factor);
% s = resample(phDiff, p, q);

% %%
% % s = sFMf0_80000;
% % envelope = s;
% % % envelope = hilbert(s) .* exp(1i * (-2*pi*f0) .* (0 : 1/fsHz : T-1/fsHz));
% a = abs(env);
% aMean = mean(a);
% aNorm = a / aMean;
% aCentNorm = aNorm - 1;
% % 
% % % Nfft = 2 ^ nextpow2(length(aCentNorm));
% % Nfft = 1024;
% % % aAbsFft = abs(fft(aCentNorm, Nfft)); % Check this ( /N or /(N/2) or nothing)
% % % gammaMax = max(aAbsFft.^2) / (Nfft/2);
% % aAbsFft = abs(fft(aCentNorm, Nfft)) / (Nfft/2); % Check this ( /N or /(N/2) or nothing)
% % gammaMax = max(aAbsFft.^2);
% % % gMax = 128754409;
% % % gMax / 16384^2;
% % 
% % aThreshold = 1;
% % phi = angle(envelope*exp(1i*(-pi/2)));
% % phiNW = phi(aNorm > aThreshold);
% % % phiUW = unwrap(phiNW);
% % phiUW = phiNW;
% % phiCent = phiUW - mean(phiUW);
% % phiN = phiCent;
% % 
% % C = length(phiN);
% % sumPhiNL2 = sum(phiN .^2);
% % sumAbsPhiNL = sum(abs(phiN));
% % sumPhiNL = sum(phiN);
% % sigmaAP = sqrt((sumPhiNL2/C) - ((sumAbsPhiNL/C) ^ 2)) % check when it is possible to become complex
% % sigmaDP = sqrt((sumPhiNL2/C) - ((sumPhiNL/C) ^ 2))
% % 
% % figure(5);
% % subplot(2,2,1); plot(angle(envelope)); grid on;
% % subplot(2,2,2); plot(phi); grid on;
% % subplot(2,2,3); plot(phiUW); grid on;
% % subplot(2,2,4); plot(phiCent); grid on;
% % 
% % figure(2);
% % subplot(2,2,1); plot(a); grid on;
% % subplot(2,2,2); plot(aNorm); grid on;
% % subplot(2,2,3); plot(aCentNorm); grid on;
% % subplot(2,2,4); plot(16384*fftshift(aAbsFft)); grid on;
% 
% 
% aThreshold = thresholds.ampl;
% phDiff = angle(env(2 : end) .* conj(env(1 : end-1)));
% aNorm = aNorm(1:end-1);
% phDiffNonWeak = phDiff(aNorm > aThreshold);
% phNWCent = phDiffNonWeak - mean(phDiffNonWeak);
% C = length(phNWCent);
% sumPh2 = sum(phNWCent .^ 2);
% sumAbsPh = sum(abs(phNWCent));
% sigmaAF = sqrt(sumPh2 / C - (sumAbsPh / C)^2)
% 
% figure(7);
% subplot(2,1,1); plot(phDiff); grid on;
% subplot(2,1,2); plot(abs(phDiff)); grid on;
% 
% figure(8);
% subplot(2,1,1); plot(phDiffNonWeak); grid on;
% subplot(2,1,2); plot(abs(phDiffNonWeak)); grid on;
% 
% figure(9);
% subplot(2,1,1); plot(phNWCent); grid on;
% subplot(2,1,2); plot(abs(phNWCent)); grid on;
