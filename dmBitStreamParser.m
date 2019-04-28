close all;
clear;

preamSym(1,:) = [3 -3  3 3  3 3  -3 -3  -3 3  3 -3  -3 3  -3 -3  3 -3  3 3  -3 -3  3 -3 ];
preamSym(2,:) = [3 -3  -3 -3  3 -3  -3 3  3 3  -3 3  -3 3  3 3  3 -3  -3 3  -3 -3  -3 3 ];
preamSym(3,:) = [3 -3  3 -3  -3 3  3 3  3 3  -3 -3  3 -3  -3 3  -3 -3  -3 3  3 -3  3 -3 ];

stream = ReadData('streamBitDmr.bin', 'int8')';
% figure(1);
% plot(stream);
% grid on;

bps = 4;
lenPackBit = 576;
lenPackSym = lenPackBit / 2;
packetsNum = floor(length(stream) / lenPackSym) / bps;
stream = stream(1 : lenPackSym * packetsNum * bps);

symbols = [stream(1:bps:end); stream(2:bps:end); stream(3:bps:end); stream(4:bps:end)];
bits = [symbols2bits(symbols(1,:)); symbols2bits(symbols(2,:)); ...
    symbols2bits(symbols(3,:)); symbols2bits(symbols(4,:))];

packsBit1 = reshape(bits(1,:)', [lenPackBit, packetsNum])';
packsBit2 = reshape(bits(2,:)', [lenPackBit, packetsNum])';
packsBit3 = reshape(bits(3,:)', [lenPackBit, packetsNum])';
packsBit4 = reshape(bits(4,:)', [lenPackBit, packetsNum])';

figure(2);
subplot(2,2,1); imagesc(packsBit1); grid on;
subplot(2,2,2); imagesc(packsBit2); grid on;
subplot(2,2,3); imagesc(packsBit3); grid on;
subplot(2,2,4); imagesc(packsBit4); grid on;

figure(3);
subplot(3,1,1); plot(xcorr(symbols(1,:), preamSym(1,:))); grid on;
subplot(3,1,2); plot(xcorr(symbols(1,:), preamSym(2,:))); grid on;
subplot(3,1,3); plot(xcorr(symbols(1,:), preamSym(3,:))); grid on;

figure(4);
subplot(3,1,1); plot(xcorr(symbols(2,:), preamSym(1,:))); grid on;
subplot(3,1,2); plot(xcorr(symbols(2,:), preamSym(2,:))); grid on;
subplot(3,1,3); plot(xcorr(symbols(2,:), preamSym(3,:))); grid on;

figure(5);
subplot(3,1,1); plot(xcorr(symbols(3,:), preamSym(1,:))); grid on;
subplot(3,1,2); plot(xcorr(symbols(3,:), preamSym(2,:))); grid on;
subplot(3,1,3); plot(xcorr(symbols(3,:), preamSym(3,:))); grid on;

figure(6);
subplot(3,1,1); plot(xcorr(symbols(4,:), preamSym(1,:))); grid on;
subplot(3,1,2); plot(xcorr(symbols(4,:), preamSym(2,:))); grid on;
subplot(3,1,3); plot(xcorr(symbols(4,:), preamSym(3,:))); grid on;

% % len = size(symbols, 2);
% len = 20;
% first = 1;
% 
% figure(2);
% subplot(4,1,1); plot(symbols(1, first:len)); grid on;
% subplot(4,1,2); plot(symbols(2, first:len)); grid on;
% subplot(4,1,3); plot(symbols(3, first:len)); grid on;
% subplot(4,1,4); plot(symbols(4, first:len)); grid on;
% 
% figure(3);
% subplot(4,1,1); plot(bits(1, first:2*len)); grid on;
% subplot(4,1,2); plot(bits(2, first:2*len)); grid on;
% subplot(4,1,3); plot(bits(3, first:2*len)); grid on;
% subplot(4,1,4); plot(bits(4, first:2*len)); grid on;