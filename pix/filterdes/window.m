bar=bartlett(100);
blk=blackman(100);
ham=hamming(100);
han=hanning(100);
plot(bar,'-')
hold on
plot(blk,'--')
plot(ham,':')
plot(han,'-.')
hold off
