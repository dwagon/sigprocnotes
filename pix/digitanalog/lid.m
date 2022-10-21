% Linear Interpolation with delay


t=.5;
f=[-2*pi:0.1:2*pi];
h=t.*((sin(pi.*f*t)./(pi.*f*t)).^2).*(exp(-j*2*pi.*f*t));
a=abs(h);
b=angle(h);
clg
subplot(211);
plot(f,a)
plot(f,b)
a=a';
b=b';
save lidfreq.data a /ascii
save lidphase.data b /ascii
