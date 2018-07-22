clear all;
a=0:(1/8000):0.5;
f40=261.626*a;
x40=sin(2*pi*f.*a);
subplot(2,1,1), plot(a,f)
subplot(2,1,2), plot(a,x_40)
audiowrite('sound1.wav',x_40,8000);
sound(x_40,8000);