%theta=-pi*f*t+atan(2*pi*f*t);
%h=t*sqrt(1+4*pi*f^2*t^2)*((sin(pi*f*t)/(pi*f*t))^2)*exp(j*theta);

t=.5;
f=[-2*pi:0.1:2*pi];
theta=-pi*f*t+atan(2*pi*f*t);
h3=exp(j*theta);
h2=(sin(pi*f*t)./(pi*f*t)).^2;
h1=sqrt(1+4*pi*f.^2*t.^2);
h=t.*h1.*h2.*h3;
a=abs(h);
b=angle(h);
subplot(211);
plot(f,a)
plot(f,b)
a=a';
b=b';
save fohfreq.data a /ascii
save fohphase.data b /ascii
