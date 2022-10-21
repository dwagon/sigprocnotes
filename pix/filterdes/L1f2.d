clear;
%format long;
%plot([boxcar(129) hamming(129) hanning(129) bartlett(129) blackman(129) chebwin(128,60)]), 
%title('Commonly used windows with length M=129'),grid
%xlabel('n'),
%ylabel('Amplitude of w[n]'), pause
%%meta windowcom1

%n=1:1:129;
%plot(n,kaiser(129,0),n,kaiser(129,1),n,kaiser(129,3),n,kaiser(129,6),n,kaiser(129,10),n,kaiser(129,20)),grid
%title('Kaiser windows with length M=129, beta=0,1,3,6,10 and 20.'),
%xlabel('n'),
%ylabel('Amplitude of w[n]'), pause
%%meta windowcom3

%a length 33 rectangular window.
winr=boxcar(33);
n=0:1:32;
winr1=winr';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winr1, [0 1 0]);
clg,
plot(n, winr1, 'o', kk, yy, '--'),
title('The rectangular window with length M=33.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude'), pause
%meta windowcom3
winr(1024)=0;
magr=abs(fft(winr));
scaler=max(magr);
magrn=magr/scaler;
windr=20*log10(magrn); %calculate in dB
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windr(1:512))
title('Magnitude frequency spectrum of the rectangular window with length M=33.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
%meta windowcom3

%a length 33 Hamming window.
winh=hamming(33);
n=0:1:32;
winh1=winh';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winh1, [0 1 0]);
clg,
plot(n, winh1, 'o', kk, yy, '--'),
title('The Hamming window with length M=33.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude'), pause
%meta windowcom3
winh(1024)=0;
magh=abs(fft(winh));
scaler=max(magh);
maghn=magh/scaler;
windh=20*log10(maghn); %calculate in dB
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windh(1:512))
title('Magnitude frequency spectrum of Hamming window with length M=33.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
%meta windowcom3

%a length 33 Blackman window.
winm=blackman(33);
n=0:1:32;
winm1=winm';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winm1, [0 1 0]);
clg,
plot(n, winm1, 'o', kk, yy, '--'),
title('The Blackman window with length M=33.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude'), pause
%meta windowcom3
winm(1024)=0;
magm=abs(fft(winm));
scaler=max(magm);
magmn=magm/scaler;
windm=20*log10(magmn); %calculate in dB
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windm(1:512))
title('Magnitude frequency spectrum of Blackman window with length M=33.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
%meta windowcom3

%a length 33 Hanning window.
winha=hanning(33);
n=0:1:32;
winha1=winha';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winha1, [0 1 0]);
clg,
plot(n, winha1, 'o', kk, yy, '--'),
title('The Hanning window with length M=33.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude'), pause
%meta windowcom3
winha(1024)=0;
magha=abs(fft(winha));
scaler=max(magha);
maghna=magha/scaler;
windha=20*log10(maghna); %calculate in dB

Warning: Log of zero
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windha(1:512))
title('Magnitude frequency spectrum of Hanning window with length M=33.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
%meta windowcom3

%a length 33 bartlett window.
winb=bartlett(33);
n=0:1:32;
winb1=winb';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winb1, [0 1 0]);
clg,
plot(n, winb1, 'o', kk, yy, '--'),
title('The Bartlett window with length M=33.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude'), pause
%meta windowcom3
winb(1024)=0;
magb=abs(fft(winb));
scaler=max(magb);
magbn=magb/scaler;
windha=20*log10(magbn); %calculate in dB

Warning: Log of zero
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windha(1:512))

Warning: Infinity or NaN found in data and not shown on plot.

title('Magnitude frequency spectrum of Bartlett window with length M=33.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
%meta windowcom3

%a length 33 Chebyshev window with r equal to 60 dB.
win=chebwin(32,60);
Chebwin: N must be odd - N is being increased by 1.
win1=win';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(win1, [0 1 0]);
clg,
plot(n, win1, 'o', kk, yy, '--'),
title('The Chebyshev window with length M=33 and r=60dB.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude'), pause
%meta windowcom3
win(1024)=0;
magc=abs(fft(win));
scaler=max(magc);
maghc=magc/scaler;
windb=20*log10(maghc); %calculate in dB
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windb(1:512))
title('Magnitude frequency spectrum of the Chebyshev window with length M=33 and r=60dB.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
%meta windowcom3

%a length 64 Hamming window.
winh2=hamming(64);
n=0:1:63;
winh21=winh2';
%k = 1:64;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winh21, [0 1 0]);
clg,
plot(n, winh21, 'o', kk, yy, '--'),
title('The Hamming window with length M=64.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude'), pause
%meta windowcom3
winh2(1024)=0;
magh2=abs(fft(winh2));
scaler=max(magh2);
maghn2=magh2/scaler;
windh2=20*log10(maghn2); %calculate in dB
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windh2(1:512))
title('Magnitude frequency spectrum of Hamming window with length M=64.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
%meta windowcom3

%fftlength = 1024;
%halfl=fftlength/2;
%%Y=fft(y,n);
%%Y=fft(y,512);
%Y=fft(y,fftlength);
%Pyy=Y.*conj(Y)/n;
%%ff=Fs*(0:(n-1)/2)/n
%%ff=Fs*(0:255)/512
%ff=Fs*(0:halfl-1)/fftlength;
%plot(ff,Pyy(1:halfl)),grid
%title('Power Spectrum, Fs=4*Fn.'),
%xlabel('Hz'),
%ylabel('Amplitude of Power Spectrum'), pause
%%plot(kk, yy, '-.', k, y, '-'), pause
%%meta samplecom

%for i=1:n;
%y(i) = cos(omega1*(i-1) - (pi/3)) + cos(omega2*(i-1)) + cos(omega3*(i-1));
%end;
%k = 0:n-1;
%yy = kron(y, [0 1 0]);
%kk = kron(k, [ 1 1 1]) ;

%tp=t*(1/T); %This is time scaling (or tp=t*Fs).
%clg,
%plot(tp,y1,'-',kk, yy, '-.', k, y, 'o'), %pause
%title('Sampling a continuous-time signal, Fs=4*Fn.'),
%xlabel('Number of samples, t*(1/T)'),
%ylabel('Amplitude'), pause

>> exit

 296858 flop(s).

clear;
%format long;
%plot([boxcar(129) hamming(129) hanning(129) bartlett(129) blackman(129) chebwin(128,60)]), 
%title('Commonly used windows with length M=129'),grid
%xlabel('n'),
%ylabel('Amplitude of w[n]'), pause
%%meta windowcom1

%n=1:1:129;
%plot(n,kaiser(129,0),n,kaiser(129,1),n,kaiser(129,3),n,kaiser(129,6),n,kaiser(129,10),n,kaiser(129,20)),grid
%title('Kaiser windows with length M=129, beta=0,1,3,6,10 and 20.'),
%xlabel('n'),
%ylabel('Amplitude of w[n]'), pause
%%meta windowcom3

%a length 33 rectangular window.
winr=boxcar(33);
n=0:1:32;
winr1=winr';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winr1, [0 1 0]);
clg,
subplot(211)
plot(n, winr1, 'o', kk, yy, '--'),
title('The rectangular window with length M=33.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude'), pause
%meta windowcom3
winr(1024)=0;
magr=abs(fft(winr));
scaler=max(magr);
magrn=magr/scaler;
windr=20*log10(magrn); %calculate in dB
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windr(1:512))
title('Magnitude frequency spectrum of the rectangular window with length M=33.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
%meta windowcom3

%a length 33 Hamming window.
winh=hamming(33);
n=0:1:32;
winh1=winh';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winh1, [0 1 0]);
clg,
subplot(211)
plot(n, winh1, 'o', kk, yy, '--'),
title('The Hamming window with length M=33.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude'), pause
%meta windowcom3
winh(1024)=0;
magh=abs(fft(winh));
scaler=max(magh);
maghn=magh/scaler;
windh=20*log10(maghn); %calculate in dB
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windh(1:512))
title('Magnitude frequency spectrum of Hamming window with length M=33.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
%meta windowcom3

%a length 33 Blackman window.
winm=blackman(33);
n=0:1:32;
winm1=winm';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winm1, [0 1 0]);
clg,
subplot(211)
plot(n, winm1, 'o', kk, yy, '--'),
title('The Blackman window with length M=33.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude'), pause
%meta windowcom3
winm(1024)=0;
magm=abs(fft(winm));
scaler=max(magm);
magmn=magm/scaler;
windm=20*log10(magmn); %calculate in dB
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windm(1:512))
title('Magnitude frequency spectrum of Blackman window with length M=33.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
%meta windowcom3

%a length 33 Hanning window.
winha=hanning(33);
n=0:1:32;
winha1=winha';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winha1, [0 1 0]);
clg,
subplot(211)
plot(n, winha1, 'o', kk, yy, '--'),
title('The Hanning window with length M=33.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude'), pause
%meta windowcom3
winha(1024)=0;
magha=abs(fft(winha));
scaler=max(magha);
maghna=magha/scaler;
windha=20*log10(maghna); %calculate in dB

Warning: Log of zero
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windha(1:512))
title('Magnitude frequency spectrum of Hanning window with length M=33.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
%meta windowcom3

%a length 33 bartlett window.
winb=bartlett(33);
n=0:1:32;
winb1=winb';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winb1, [0 1 0]);
clg,
subplot(211)
plot(n, winb1, 'o', kk, yy, '--'),
title('The Bartlett window with length M=33.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude'), pause
%meta windowcom3
winb(1024)=0;
magb=abs(fft(winb));
scaler=max(magb);
magbn=magb/scaler;
windha=20*log10(magbn); %calculate in dB

Warning: Log of zero
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windha(1:512))

Warning: Infinity or NaN found in data and not shown on plot.

title('Magnitude frequency spectrum of Bartlett window with length M=33.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
%meta windowcom3

%a length 33 Chebyshev window with r equal to 60 dB.
win=chebwin(32,60);
Chebwin: N must be odd - N is being increased by 1.
win1=win';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(win1, [0 1 0]);
clg,
subplot(211)
plot(n, win1, 'o', kk, yy, '--'),
title('The Chebyshev window with length M=33 and r=60dB.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude'), pause
%meta windowcom3
win(1024)=0;
magc=abs(fft(win));
scaler=max(magc);
maghc=magc/scaler;
windb=20*log10(maghc); %calculate in dB
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windb(1:512))
title('Magnitude frequency spectrum of the Chebyshev window with length M=33 and r=60dB.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
%meta windowcom3

%a length 64 Hamming window.
winh2=hamming(64);
n=0:1:63;
winh21=winh2';
%k = 1:64;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winh21, [0 1 0]);
clg,
subplot(211)
plot(n, winh21, 'o', kk, yy, '--'),
title('The Hamming window with length M=64.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude'), pause
%meta windowcom3
winh2(1024)=0;
magh2=abs(fft(winh2));
scaler=max(magh2);
maghn2=magh2/scaler;
windh2=20*log10(maghn2); %calculate in dB
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windh2(1:512))
title('Magnitude frequency spectrum of Hamming window with length M=64.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
%meta windowcom3

%fftlength = 1024;
%halfl=fftlength/2;
%%Y=fft(y,n);
%%Y=fft(y,512);
%Y=fft(y,fftlength);
%Pyy=Y.*conj(Y)/n;
%%ff=Fs*(0:(n-1)/2)/n
%%ff=Fs*(0:255)/512
%ff=Fs*(0:halfl-1)/fftlength;
%plot(ff,Pyy(1:halfl)),grid
%title('Power Spectrum, Fs=4*Fn.'),
%xlabel('Hz'),
%ylabel('Amplitude of Power Spectrum'), pause
%%plot(kk, yy, '-.', k, y, '-'), pause
%%meta samplecom

%for i=1:n;
%y(i) = cos(omega1*(i-1) - (pi/3)) + cos(omega2*(i-1)) + cos(omega3*(i-1));
%end;
%k = 0:n-1;
%yy = kron(y, [0 1 0]);
%kk = kron(k, [ 1 1 1]) ;

%tp=t*(1/T); %This is time scaling (or tp=t*Fs).
%clg,
%plot(tp,y1,'-',kk, yy, '-.', k, y, 'o'), %pause
%title('Sampling a continuous-time signal, Fs=4*Fn.'),
%xlabel('Number of samples, t*(1/T)'),
%ylabel('Amplitude'), pause

>> window3
%This is a sampling plot and it not quite right.
diary L1f2.d 
echo on
clear;
%format long;
%plot([boxcar(129) hamming(129) hanning(129) bartlett(129) blackman(129) chebwin(128,60)]), 
%title('Commonly used windows with length M=129'),grid
%xlabel('n'),
%ylabel('Amplitude of w[n]'), pause
%%meta windowcom1

%n=1:1:129;
%plot(n,kaiser(129,0),n,kaiser(129,1),n,kaiser(129,3),n,kaiser(129,6),n,kaiser(129,10),n,kaiser(129,20)),grid
%title('Kaiser windows with length M=129, beta=0,1,3,6,10 and 20.'),
%xlabel('n'),
%ylabel('Amplitude of w[n]'), pause
%%meta windowcom3

%a length 33 rectangular window.
winr=boxcar(33);
n=0:1:32;
winr1=winr';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winr1, [0 1 0]);
clg,
subplot(211)
plot(n, winr1, 'o', kk, yy, '--'),
title('The rectangular window with length M=33.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude')
%meta windowcom3
winr(1024)=0;
magr=abs(fft(winr));
scaler=max(magr);
magrn=magr/scaler;
windr=20*log10(magrn); %calculate in dB
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windr(1:512))
title('Magnitude frequency spectrum of the rectangular window with length M=33.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
%meta windowcom3

%a length 33 Hamming window.
winh=hamming(33);
n=0:1:32;
winh1=winh';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winh1, [0 1 0]);
clg,
subplot(211)
plot(n, winh1, 'o', kk, yy, '--'),
title('The Hamming window with length M=33.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude')
%meta windowcom3
winh(1024)=0;
magh=abs(fft(winh));
scaler=max(magh);
maghn=magh/scaler;
windh=20*log10(maghn); %calculate in dB
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windh(1:512))
title('Magnitude frequency spectrum of Hamming window with length M=33.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
%meta windowcom3

%a length 33 Blackman window.
winm=blackman(33);
n=0:1:32;
winm1=winm';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winm1, [0 1 0]);
clg,
subplot(211)
plot(n, winm1, 'o', kk, yy, '--'),
title('The Blackman window with length M=33.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude')
%meta windowcom3
winm(1024)=0;
magm=abs(fft(winm));
scaler=max(magm);
magmn=magm/scaler;
windm=20*log10(magmn); %calculate in dB
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windm(1:512))
title('Magnitude frequency spectrum of Blackman window with length M=33.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
%meta windowcom3

%a length 33 Hanning window.
winha=hanning(33);
n=0:1:32;
winha1=winha';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winha1, [0 1 0]);
clg,
subplot(211)
plot(n, winha1, 'o', kk, yy, '--'),
title('The Hanning window with length M=33.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude')
%meta windowcom3
winha(1024)=0;
magha=abs(fft(winha));
scaler=max(magha);
maghna=magha/scaler;
windha=20*log10(maghna); %calculate in dB

Warning: Log of zero
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windha(1:512))
title('Magnitude frequency spectrum of Hanning window with length M=33.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
%meta windowcom3

%a length 33 bartlett window.
winb=bartlett(33);
n=0:1:32;
winb1=winb';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winb1, [0 1 0]);
clg,
subplot(211)
plot(n, winb1, 'o', kk, yy, '--'),
title('The Bartlett window with length M=33.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude')
%meta windowcom3
winb(1024)=0;
magb=abs(fft(winb));
scaler=max(magb);
magbn=magb/scaler;
windha=20*log10(magbn); %calculate in dB

Warning: Log of zero
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windha(1:512))

Warning: Infinity or NaN found in data and not shown on plot.

title('Magnitude frequency spectrum of Bartlett window with length M=33.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
%meta windowcom3

%a length 33 Chebyshev window with r equal to 60 dB.
win=chebwin(32,60);
Chebwin: N must be odd - N is being increased by 1.
win1=win';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(win1, [0 1 0]);
clg,
subplot(211)
plot(n, win1, 'o', kk, yy, '--'),
title('The Chebyshev window with length M=33 and r=60dB.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude')
%meta windowcom3
win(1024)=0;
magc=abs(fft(win));
scaler=max(magc);
maghc=magc/scaler;
windb=20*log10(maghc); %calculate in dB
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windb(1:512))
title('Magnitude frequency spectrum of the Chebyshev window with length M=33 and r=60dB.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
%meta windowcom3

%a length 64 Hamming window.
winh2=hamming(64);
n=0:1:63;
winh21=winh2';
%k = 1:64;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winh21, [0 1 0]);
clg,
subplot(211)
plot(n, winh21, 'o', kk, yy, '--'),
title('The Hamming window with length M=64.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude')
%meta windowcom3
winh2(1024)=0;
magh2=abs(fft(winh2));
scaler=max(magh2);
maghn2=magh2/scaler;
windh2=20*log10(maghn2); %calculate in dB
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windh2(1:512))
title('Magnitude frequency spectrum of Hamming window with length M=64.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
%meta windowcom3

%fftlength = 1024;
%halfl=fftlength/2;
%%Y=fft(y,n);
%%Y=fft(y,512);
%Y=fft(y,fftlength);
%Pyy=Y.*conj(Y)/n;
%%ff=Fs*(0:(n-1)/2)/n
%%ff=Fs*(0:255)/512
%ff=Fs*(0:halfl-1)/fftlength;
%plot(ff,Pyy(1:halfl)),grid
%title('Power Spectrum, Fs=4*Fn.'),
%xlabel('Hz'),
%ylabel('Amplitude of Power Spectrum'), pause
%%plot(kk, yy, '-.', k, y, '-'), pause
%%meta samplecom

%for i=1:n;
%y(i) = cos(omega1*(i-1) - (pi/3)) + cos(omega2*(i-1)) + cos(omega3*(i-1));
%end;
%k = 0:n-1;
%yy = kron(y, [0 1 0]);
%kk = kron(k, [ 1 1 1]) ;

%tp=t*(1/T); %This is time scaling (or tp=t*Fs).
%clg,
%plot(tp,y1,'-',kk, yy, '-.', k, y, 'o'), %pause
%title('Sampling a continuous-time signal, Fs=4*Fn.'),
%xlabel('Number of samples, t*(1/T)'),
%ylabel('Amplitude'), pause

>> 
>> window3
%This is a sampling plot and it not quite right.
diary L1f2.d 
echo on
clear;
%format long;
%plot([boxcar(129) hamming(129) hanning(129) bartlett(129) blackman(129) chebwin(128,60)]), 
%title('Commonly used windows with length M=129'),grid
%xlabel('n'),
%ylabel('Amplitude of w[n]'), pause
%%meta windowcom1

%n=1:1:129;
%plot(n,kaiser(129,0),n,kaiser(129,1),n,kaiser(129,3),n,kaiser(129,6),n,kaiser(129,10),n,kaiser(129,20)),grid
%title('Kaiser windows with length M=129, beta=0,1,3,6,10 and 20.'),
%xlabel('n'),
%ylabel('Amplitude of w[n]'), pause
%%meta windowcom3

%a length 33 rectangular window.
winr=boxcar(33);
n=0:1:32;
winr1=winr';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winr1, [0 1 0]);
clg,
subplot(211)
plot(n, winr1, 'o', kk, yy, '--'),
title('The rectangular window with length M=33.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude')
winr(1024)=0;
magr=abs(fft(winr));
scaler=max(magr);
magrn=magr/scaler;
windr=20*log10(magrn); %calculate in dB
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windr(1:512))
title('Magnitude frequency spectrum of the rectangular window with length M=33.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
meta rectwind

%a length 33 Hamming window.
winh=hamming(33);
n=0:1:32;
winh1=winh';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winh1, [0 1 0]);
clg,
subplot(211)
plot(n, winh1, 'o', kk, yy, '--'),
title('The Hamming window with length M=33.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude')
winh(1024)=0;
magh=abs(fft(winh));
scaler=max(magh);
maghn=magh/scaler;
windh=20*log10(maghn); %calculate in dB
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windh(1:512))
title('Magnitude frequency spectrum of Hamming window with length M=33.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
meta hammwind

%a length 33 Blackman window.
winm=blackman(33);
n=0:1:32;
winm1=winm';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winm1, [0 1 0]);
clg,
subplot(211)
plot(n, winm1, 'o', kk, yy, '--'),
title('The Blackman window with length M=33.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude')
winm(1024)=0;
magm=abs(fft(winm));
scaler=max(magm);
magmn=magm/scaler;
windm=20*log10(magmn); %calculate in dB
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windm(1:512))
title('Magnitude frequency spectrum of Blackman window with length M=33.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
meta blackwind

%a length 33 Hanning window.
winha=hanning(33);
n=0:1:32;
winha1=winha';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winha1, [0 1 0]);
clg,
subplot(211)
plot(n, winha1, 'o', kk, yy, '--'),
title('The Hanning window with length M=33.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude')
winha(1024)=0;
magha=abs(fft(winha));
scaler=max(magha);
maghna=magha/scaler;
windha=20*log10(maghna); %calculate in dB

Warning: Log of zero
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windha(1:512))
title('Magnitude frequency spectrum of Hanning window with length M=33.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
meta hannwind

%a length 33 bartlett window.
winb=bartlett(33);
n=0:1:32;
winb1=winb';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winb1, [0 1 0]);
clg,
subplot(211)
plot(n, winb1, 'o', kk, yy, '--'),
title('The Bartlett window with length M=33.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude')
winb(1024)=0;
magb=abs(fft(winb));
scaler=max(magb);
magbn=magb/scaler;
windha=20*log10(magbn); %calculate in dB

Warning: Log of zero
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windha(1:512))

Warning: Infinity or NaN found in data and not shown on plot.

title('Magnitude frequency spectrum of Bartlett window with length M=33.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
meta bartwind

%a length 33 Chebyshev window with r equal to 60 dB.
win=chebwin(32,60);
Chebwin: N must be odd - N is being increased by 1.
win1=win';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(win1, [0 1 0]);
clg,
subplot(211)
plot(n, win1, 'o', kk, yy, '--'),
title('The Chebyshev window with length M=33 and r=60dB.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude')
win(1024)=0;
magc=abs(fft(win));
scaler=max(magc);
maghc=magc/scaler;
windb=20*log10(maghc); %calculate in dB
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windb(1:512))
title('Magnitude frequency spectrum of the Chebyshev window with length M=33 and r=60dB.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
meta chebywind

%a length 64 Hamming window.
winh2=hamming(64);
n=0:1:63;
winh21=winh2';
%k = 1:64;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winh21, [0 1 0]);
clg,
subplot(211)
plot(n, winh21, 'o', kk, yy, '--'),
title('The Hamming window with length M=64.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude')
winh2(1024)=0;
magh2=abs(fft(winh2));
scaler=max(magh2);
maghn2=magh2/scaler;
windh2=20*log10(maghn2); %calculate in dB
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windh2(1:512))
title('Magnitude frequency spectrum of Hamming window with length M=64.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
meta hamm64wind

%fftlength = 1024;
%halfl=fftlength/2;
%%Y=fft(y,n);
%%Y=fft(y,512);
%Y=fft(y,fftlength);
%Pyy=Y.*conj(Y)/n;
%%ff=Fs*(0:(n-1)/2)/n
%%ff=Fs*(0:255)/512
%ff=Fs*(0:halfl-1)/fftlength;
%plot(ff,Pyy(1:halfl)),grid
%title('Power Spectrum, Fs=4*Fn.'),
%xlabel('Hz'),
%ylabel('Amplitude of Power Spectrum'), pause
%%plot(kk, yy, '-.', k, y, '-'), pause
%%meta samplecom

%for i=1:n;
%y(i) = cos(omega1*(i-1) - (pi/3)) + cos(omega2*(i-1)) + cos(omega3*(i-1));
%end;
%k = 0:n-1;
%yy = kron(y, [0 1 0]);
%kk = kron(k, [ 1 1 1]) ;

%tp=t*(1/T); %This is time scaling (or tp=t*Fs).
%clg,
%plot(tp,y1,'-',kk, yy, '-.', k, y, 'o'), %pause
%title('Sampling a continuous-time signal, Fs=4*Fn.'),
%xlabel('Number of samples, t*(1/T)'),
%ylabel('Amplitude'), pause

>> meta
>> exit

 890574 flop(s).

clear;
%format long;
%plot([boxcar(129) hamming(129) hanning(129) bartlett(129) blackman(129) chebwin(128,60)]), 
%title('Commonly used windows with length M=129'),grid
%xlabel('n'),
%ylabel('Amplitude of w[n]'), pause
%%meta windowcom1

%n=1:1:129;
%plot(n,kaiser(129,0),n,kaiser(129,1),n,kaiser(129,3),n,kaiser(129,6),n,kaiser(129,10),n,kaiser(129,20)),grid
%title('Kaiser windows with length M=129, beta=0,1,3,6,10 and 20.'),
%xlabel('n'),
%ylabel('Amplitude of w[n]'), pause
%%meta windowcom3

%a length 33 rectangular window.
winr=boxcar(33);
n=0:1:32;
winr1=winr';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winr1, [0 1 0]);
clg,
subplot(211)
plot(n, winr1, 'o', kk, yy, '--'),
title('The rectangular window with length M=33.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude')
winr(1024)=0;
magr=abs(fft(winr));
scaler=max(magr);
magrn=magr/scaler;
windr=20*log10(magrn); %calculate in dB
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windr(1:512))
title('Magnitude frequency spectrum of the rectangular window with length M=33.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
%meta rectwind

%a length 33 Hamming window.
winh=hamming(33);
n=0:1:32;
winh1=winh';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winh1, [0 1 0]);
clg,
subplot(211)
plot(n, winh1, 'o', kk, yy, '--'),
title('The Hamming window with length M=33.'),
xlabel('Length of the window, M.'),
>> kaiser
echo off
??? plot(n,kaiser(
                 |
Improper function reference. A "," or ")" is expected.

Error in ==> kaiser.m
On line 5 ==> plot(n,kaiser(129,0),n,kaiser(129,1),n,kaiser(129,3),n,kaiser(129,6),n,kaiser(129,10),n,kaiser(129,20)),grid

>> help kaiser

 meta windowcom3

>> kais
??? plot(n,kaiser(
                 |
Improper function reference. A "," or ")" is expected.

Error in ==> kais.m
On line 5 ==> plot(n,kaiser(129,0),n,kaiser(129,1),n,kaiser(129,3),n,kaiser(129,6),n,kaiser(129,10),n,kaiser(129,20)),grid

>> help kaiser

 KAISER	 KAISER(N,beta) returns the BETA-valued N-point Kaiser window.

>> kaiser
???  Input argument is undefined.
Symbol in question ==> nn

Error in ==> /matlab/matlab/toolbox/signal/kaiser.m
On line 7  ==> nw = round(nn);

>> winkais.m
??? winkais.
           |
Missing operator, comma, or semi-colon.

>> clear
>> winkais
???  Undefined function or variable.
Symbol in question ==> winr1

Error in ==> winkais.m
On line 11 ==> plot(n, winr1, 'o', kk, yy, '--'),

>> winkais
???  Undefined function or variable.
Symbol in question ==> winr1

Error in ==> winkais.m
On line 12 ==> plot(n, winr1, 'o', kk, yy, '--'),

>> window3
clear;
%format long;
%plot([boxcar(129) hamming(129) hanning(129) bartlett(129) blackman(129) chebwin(128,60)]), 
%title('Commonly used windows with length M=129'),grid
%xlabel('n'),
%ylabel('Amplitude of w[n]'), pause
%%meta windowcom1

n=1:1:129;
plot(n,kaiser(129,0),n,kaiser(129,1),n,kaiser(129,3),n,kaiser(129,6),n,kaiser(129,10),n,kaiser(129,20)),grid
title('Kaiser windows with length M=129, beta=0,1,3,6,10 and 20.'),
xlabel('n'),
ylabel('Amplitude of w[n]'), pause
clg,
subplot(211)
plot(n, winr1, 'o', kk, yy, '--'),
???  Undefined function or variable.
Symbol in question ==> winr1

>> window3
%This is a sampling plot and it not quite right.
diary L1f2.d 
echo on
clear;
%format long;
%plot([boxcar(129) hamming(129) hanning(129) bartlett(129) blackman(129) chebwin(128,60)]), 
%title('Commonly used windows with length M=129'),grid
%xlabel('n'),
%ylabel('Amplitude of w[n]'), pause
%%meta windowcom1

%n=1:1:129;
%plot(n,kaiser(129,0),n,kaiser(129,1),n,kaiser(129,3),n,kaiser(129,6),n,kaiser(129,10),n,kaiser(129,20)),grid
%title('Kaiser windows with length M=129, beta=0,1,3,6,10 and 20.'),
%xlabel('n'),
%ylabel('Amplitude of w[n]'), pause
%clg,
%subplot(211)
%plot(n, winr1, 'o', kk, yy, '--'),
%title('The rectangular window with length M=33.'),
%xlabel('Length of the window, M.'),
%ylabel('Magnitude')
%winr(1024)=0;
%magr=abs(fft(winr));
%scaler=max(magr);
%magrn=magr/scaler;
%windr=20*log10(magrn); %calculate in dB
%n1=1:1:512;
%fr=n1/(2*512);
%plot(fr,windr(1:512))
%title('Magnitude frequency spectrum of the rectangular window with length M=33.'),
%xlabel('Relative frequency f in cycles per sample'),
%ylabel('Magnitude spectrum in dB'), pause
%%meta windowcom3

%a length 33 rectangular window.
winr=boxcar(33);
n=0:1:32;
winr1=winr';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winr1, [0 1 0]);
clg,
subplot(211)
plot(n, winr1, 'o', kk, yy, '--'),
title('The rectangular window with length M=33.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude')
winr(1024)=0;
magr=abs(fft(winr));
scaler=max(magr);
magrn=magr/scaler;
windr=20*log10(magrn); %calculate in dB
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windr(1:512))
title('Magnitude frequency spectrum of the rectangular window with length M=33.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause

>> help kron

 KRON	KRON(X,Y) is the Kronecker tensor product of X and Y.
 	The result is a large matrix formed by taking all possible
 	products between the elements of X and those of Y.   For
 	example, if X is 2 by 3, then KRON(X,Y) is
 
 	   [ X(1,1)*Y  X(1,2)*Y  X(1,3)*Y
 	     X(2,1)*Y  X(2,2)*Y  X(2,3)*Y ]

>> 
>> winkais
echo off
>> help kasier


Sorry, no HELP on kasier

>> help kaiser

 KAISER	 KAISER(N,beta) returns the BETA-valued N-point Kaiser window.

>> winkais
>> whos
                 Name        Size          Total     Complex

                   fr       1 by 512         512        No
                   kk       1 by 387         387        No
                 magr    1024 by 1          1024        No
                magrn    1024 by 1          1024        No
                    n       1 by 129         129        No
                   n1       1 by 512         512        No
               scaler       1 by 1             1        No
                   w1    1024 by 1          1024        No
                  w1t       1 by 129         129        No
                   w2     129 by 1           129        No
                  w2t       1 by 129         129        No
                windr    1024 by 1          1024        No
                   yy       1 by 387         387        No

Grand total is (6411 * 8) = 51288 bytes,

>> plot(n,kaiser(129,0))
>> plot(n,kaiser(129,0)')
>> plot(n,kaiser(129,1))
>> clg
>> winkais
??? ctrum in
          |
Missing operator, comma, or semi-colon.

Error in ==> winkais.m
On line 37 ==> ctrum in dB'), pause

>> winkais
???  Undefined function or variable.
Symbol in question ==> kaiswind1

Error in ==> winkais.m
On line 38 ==> kaiswind1

>> winkais
??? pectrum of
            |
Missing operator, comma, or semi-colon.

Error in ==> winkais.m
On line 38 ==> pectrum of the Kaiser window with length M=129.  Beta=0 and 1'),

>> winkais
>> help plot

PLOT	Plot vectors or matrices. PLOT(X,Y) plots vector X versus
	vector Y. If X or Y is a matrix, then the vector is plotted
	versus the rows or columns of the matrix, whichever lines
	up. PLOT(X1,Y1,X2,Y2) is another way of producing multiple
	lines on the plot. PLOT(X1,Y1,':',X2,Y2,'+') uses a
	dotted line for the first curve and the point symbol + 
	for the second curve. Other line and point types are:

	   solid   -      point  .      red       r
	   dashed  --     plus   +      green     g
	   dotted  :      star   *      blue      b
	   dashdot -.     circle o      white     w
	                  x-mark x      invisible i
	                                arbitrary c1, c15, etc.

	PLOT(Y) plots the columns of Y versus their index. PLOT(Y)
	is equivalent to PLOT(real(Y),imag(Y)) if Y is complex.
	In all other uses of PLOT, the imaginary part is ignored.
	See SEMI, LOGLOG, POLAR, GRID, SHG, CLC, CLG, TITLE, XLABEL
	YLABEL, AXIS, HOLD, MESH, CONTOUR, SUBPLOT.

>> winkais
>> flops

ans =

     1323645

>> winkais
>> flops

ans =

     1732712

>> winkais

Error in ==> winkais.m
On line 43 ==> ylabel('Magnitude spectrum in dB'), pause

>> help plot

PLOT	Plot vectors or matrices. PLOT(X,Y) plots vector X versus
	vector Y. If X or Y is a matrix, then the vector is plotted
	versus the rows or columns of the matrix, whichever lines
	up. PLOT(X1,Y1,X2,Y2) is another way of producing multiple
	lines on the plot. PLOT(X1,Y1,':',X2,Y2,'+') uses a
	dotted line for the first curve and the point symbol + 
	for the second curve. Other line and point types are:

	   solid   -      point  .      red       r
	   dashed  --     plus   +      green     g
	   dotted  :      star   *      blue      b
	   dashdot -.     circle o      white     w
	                  x-mark x      invisible i
	                                arbitrary c1, c15, etc.

	PLOT(Y) plots the columns of Y versus their index. PLOT(Y)
	is equivalent to PLOT(real(Y),imag(Y)) if Y is complex.
	In all other uses of PLOT, the imaginary part is ignored.
	See SEMI, LOGLOG, POLAR, GRID, SHG, CLC, CLG, TITLE, XLABEL
	YLABEL, AXIS, HOLD, MESH, CONTOUR, SUBPLOT.

>> winkais
>> winkais
>> clg
>> plot(fr,windr2(1:512),'-')
>> hold

Current plot held

>> plot(fr,windr(1:512),'-'),hold
>> winkais
>> winkais
>> meta
>> flops

ans =

     3482808

>> exit

 3482808 flop(s).

clear;
%format long;
%plot([boxcar(129) hamming(129) hanning(129) bartlett(129) blackman(129) chebwin(128,60)]), 
%title('Commonly used windows with length M=129'),grid
%xlabel('n'),
%ylabel('Amplitude of w[n]'), pause
%%meta windowcom1

%n=1:1:129;
%plot(n,kaiser(129,0),n,kaiser(129,1),n,kaiser(129,3),n,kaiser(129,6),n,kaiser(129,10),n,kaiser(129,20)),grid
%title('Kaiser windows with length M=129, beta=0,1,3,6,10 and 20.'),
%xlabel('n'),
%ylabel('Amplitude of w[n]'), pause
%clg,
%subplot(211)
%plot(n, winr1, 'o', kk, yy, '--'),
%title('The rectangular window with length M=33.'),
%xlabel('Length of the window, M.'),
%ylabel('Magnitude')
%winr(1024)=0;
%magr=abs(fft(winr));
%scaler=max(magr);
%magrn=magr/scaler;
%windr=20*log10(magrn); %calculate in dB
%n1=1:1:512;
%fr=n1/(2*512);
%plot(fr,windr(1:512))
%title('Magnitude frequency spectrum of the rectangular window with length M=33.'),
%xlabel('Relative frequency f in cycles per sample'),
%ylabel('Magnitude spectrum in dB'), pause
%%meta windowcom3

%a length 33 rectangular window.
winr=boxcar(33);
n=0:1:32;
winr1=winr';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winr1, [0 1 0]);
clg,
subplot(211)
plot(n, winr1, 'o', kk, yy, '--'),
title('The rectangular window with length M=33.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude')
winr(1024)=0;
magr=abs(fft(winr));
scaler=max(magr);
magrn=magr/scaler;
windr=20*log10(magrn); %calculate in dB
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windr(1:512))
title('Magnitude frequency spectrum of the rectangular window with length M=33.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
meta rectwind

%a length 33 Hamming window.
winh=hamming(33);
n=0:1:32;
winh1=winh';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winh1, [0 1 0]);
clg,
subplot(211)
plot(n, winh1, 'o', kk, yy, '--'),
title('The Hamming window with length M=33.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude')
winh(1024)=0;
magh=abs(fft(winh));
scaler=max(magh);
maghn=magh/scaler;
windh=20*log10(maghn); %calculate in dB
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windh(1:512))
title('Magnitude frequency spectrum of Hamming window with length M=33.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
meta hammwind

%a length 33 Blackman window.
winm=blackman(33);
n=0:1:32;
winm1=winm';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winm1, [0 1 0]);
clg,
subplot(211)
plot(n, winm1, 'o', kk, yy, '--'),
title('The Blackman window with length M=33.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude')
winm(1024)=0;
magm=abs(fft(winm));
scaler=max(magm);
magmn=magm/scaler;
windm=20*log10(magmn); %calculate in dB
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windm(1:512))
title('Magnitude frequency spectrum of Blackman window with length M=33.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
meta blackwind

%a length 33 Hanning window.
winha=hanning(33);
n=0:1:32;
winha1=winha';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winha1, [0 1 0]);
clg,
subplot(211)
plot(n, winha1, 'o', kk, yy, '--'),
title('The Hanning window with length M=33.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude')
winha(1024)=0;
magha=abs(fft(winha));
scaler=max(magha);
maghna=magha/scaler;
windha=20*log10(maghna); %calculate in dB

Warning: Log of zero
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windha(1:512))
title('Magnitude frequency spectrum of Hanning window with length M=33.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
meta hannwind

%a length 33 bartlett window.
winb=bartlett(33);
n=0:1:32;
winb1=winb';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winb1, [0 1 0]);
clg,
subplot(211)
plot(n, winb1, 'o', kk, yy, '--'),
title('The Bartlett window with length M=33.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude')
winb(1024)=0;
magb=abs(fft(winb));
scaler=max(magb);
magbn=magb/scaler;
windha=20*log10(magbn); %calculate in dB

Warning: Log of zero
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windha(1:512))

Warning: Infinity or NaN found in data and not shown on plot.

title('Magnitude frequency spectrum of Bartlett window with length M=33.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
meta bartwind

%a length 33 Chebyshev window with r equal to 60 dB.
win=chebwin(32,60);
Chebwin: N must be odd - N is being increased by 1.
win1=win';
%k = 1:33;
kk = kron(n, [ 1 1 1]) ;
yy = kron(win1, [0 1 0]);
clg,
subplot(211)
plot(n, win1, 'o', kk, yy, '--'),
title('The Chebyshev window with length M=33 and r=60dB.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude')
win(1024)=0;
magc=abs(fft(win));
scaler=max(magc);
maghc=magc/scaler;
windb=20*log10(maghc); %calculate in dB
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windb(1:512))
title('Magnitude frequency spectrum of the Chebyshev window with length M=33 and r=60dB.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
meta chebywind

%a length 64 Hamming window.
winh2=hamming(64);
n=0:1:63;
winh21=winh2';
%k = 1:64;
kk = kron(n, [ 1 1 1]) ;
yy = kron(winh21, [0 1 0]);
clg,
subplot(211)
plot(n, winh21, 'o', kk, yy, '--'),
title('The Hamming window with length M=64.'),
xlabel('Length of the window, M.'),
ylabel('Magnitude')
winh2(1024)=0;
magh2=abs(fft(winh2));
scaler=max(magh2);
maghn2=magh2/scaler;
windh2=20*log10(maghn2); %calculate in dB
n1=1:1:512;
fr=n1/(2*512);
plot(fr,windh2(1:512))
title('Magnitude frequency spectrum of Hamming window with length M=64.'),
xlabel('Relative frequency f in cycles per sample'),
ylabel('Magnitude spectrum in dB'), pause
meta hamm64wind

%fftlength = 1024;
%halfl=fftlength/2;
%%Y=fft(y,n);
%%Y=fft(y,512);
%Y=fft(y,fftlength);
%Pyy=Y.*conj(Y)/n;
%%ff=Fs*(0:(n-1)/2)/n
%%ff=Fs*(0:255)/512
%ff=Fs*(0:halfl-1)/fftlength;
%plot(ff,Pyy(1:halfl)),grid
%title('Power Spectrum, Fs=4*Fn.'),
%xlabel('Hz'),
%ylabel('Amplitude of Power Spectrum'), pause
%%plot(kk, yy, '-.', k, y, '-'), pause
%%meta samplecom

%for i=1:n;
%y(i) = cos(omega1*(i-1) - (pi/3)) + cos(omega2*(i-1)) + cos(omega3*(i-1));
%end;
%k = 0:n-1;
%yy = kron(y, [0 1 0]);
%kk = kron(k, [ 1 1 1]) ;

%tp=t*(1/T); %This is time scaling (or tp=t*Fs).
%clg,
%plot(tp,y1,'-',kk, yy, '-.', k, y, 'o'), %pause
%title('Sampling a continuous-time signal, Fs=4*Fn.'),
%xlabel('Number of samples, t*(1/T)'),
%ylabel('Amplitude'), pause

>> winkais
echo off
>> flops

ans =

      705925

>> exit

 705925 flop(s).

