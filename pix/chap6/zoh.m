clg;
t=1;
f=[-2*pi:0.1:2*pi];
h2=sin(pi.*f*t)./(pi.*f*t);
h3=exp(-j*pi.*f*t);
h=t.*h2.*h3;
a=abs(h);
b=angle(h);
subplot(211);
plot(f,a)
plot(f,b)
a=a';
b=b';
save zohfreq.data a /ascii
save zohphase.data b /ascii
