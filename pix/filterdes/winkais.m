echo off
!rm -f kaiswind1 kaiswind2 kaiswind3
clear;
clg;

% Kaiser with beta=0, beta=1
n=1:1:129;
w1=kaiser(129,0);
w1t=w1';
kk = kron(n,[1 1 1]);
yy=kron(w1t,[0 1 0]);
w2=kaiser(129,1);
w2t=w2';
kk2 = kron(n,[1 1 1]);
yy2=kron(w2t,[0 1 0]);
%plot(n,w1,n,w2),grid
%title('Kaiser windows with length M=129, beta=0 and 1.'),
%xlabel('n'),
%ylabel('Amplitude of w[n]'), pause
clg,
subplot(211)
plot(n, w1t, 'o', kk, yy, '--'), hold
plot(n, w2t, 'x', kk2, yy2, '--'), hold
title('The Kaiser window with length M=129. Beta=0 and 1'),
xlabel('Length of the window, M.'),
ylabel('Magnitude')
w1(1024)=0;
w2(1024)=0;
magr=abs(fft(w1));
magr2=abs(fft(w2));
scaler=max(magr);
scaler2=max(magr2);
magrn=magr/scaler;
magrn2=magr2/scaler2;
windr=20*log10(magrn); %calculate in dB
windr2=20*log10(magrn2); %calculate in dB
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windr(1:512),'-'),hold
plot(fr,windr2(1:512),'--'),hold
title('Magnitude frequency spectrum of the Kaiser window with length M=129.  Beta=0 and 1'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
meta kaiswind1

% Kaiser with beta=3, beta=6
n=1:1:129;
w1=kaiser(129,3);
w1t=w1';
kk = kron(n,[1 1 1]);
yy=kron(w1t,[0 1 0]);
w2=kaiser(129,6);
w2t=w2';
kk2 = kron(n,[1 1 1]);
yy2=kron(w2t,[0 1 0]);
%plot(n,w1,n,w2),grid
%title('Kaiser windows with length M=129, beta=3 and 6.'),
%xlabel('n'),
%ylabel('Amplitude of w[n]'), pause
clg,
subplot(211)
plot(n, w1t, 'o', kk, yy, '--'), hold
plot(n, w2t, 'x', kk2, yy2, '--'), hold
title('The Kaiser window with length M=129. Beta=3 and 6'),
xlabel('Length of the window, M.'),
ylabel('Magnitude')
w1(1024)=0;
w2(1024)=0;
magr=abs(fft(w1));
magr2=abs(fft(w2));
scaler=max(magr);
scaler2=max(magr2);
magrn=magr/scaler;
magrn2=magr2/scaler2;
windr=20*log10(magrn); %calculate in dB
windr2=20*log10(magrn2); %calculate in dB
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windr(1:512),'-'),hold
plot(fr,windr2(1:512),'--'),hold
title('Magnitude frequency spectrum of the Kaiser window with length M=129.  Beta=3 and 6'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
meta kaiswind2

% Kaiser with beta=10, beta=20
n=1:1:129;
w1=kaiser(129,10);
w1t=w1';
kk = kron(n,[1 1 1]);
yy=kron(w1t,[0 1 0]);
w2=kaiser(129,20);
w2t=w2';
kk2 = kron(n,[1 1 1]);
yy2=kron(w2t,[0 1 0]);
%plot(n,w1,n,w2),grid
%title('Kaiser windows with length M=129, beta=10 and 20.'),
%xlabel('n'),
%ylabel('Amplitude of w[n]'), pause
clg,
subplot(211)
plot(n, w1t, 'o', kk, yy, '--'), hold
plot(n, w2t, 'x', kk2, yy2, '--'), hold
title('The Kaiser window with length M=129. Beta=10 and 20'),
xlabel('Length of the window, M.'),
ylabel('Magnitude')
w1(1024)=0;
w2(1024)=0;
magr=abs(fft(w1));
magr2=abs(fft(w2));
scaler=max(magr);
scaler2=max(magr2);
magrn=magr/scaler;
magrn2=magr2/scaler2;
windr=20*log10(magrn); %calculate in dB
windr2=20*log10(magrn2); %calculate in dB
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windr2(1:512),'--'),hold
plot(fr,windr(1:512),'-'),hold
title('Magnitude frequency spectrum of the Kaiser window with length M=129.  Beta=10 and 20'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
meta kaiswind3
