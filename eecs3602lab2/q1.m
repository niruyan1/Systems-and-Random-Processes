clear all;
a=0:(1/8000):2;
f=150*a;
x=sin(2*pi*f.*a);
subplot(2,1,1), plot(a,f)
subplot(2,1,2), plot(a,x)
audiowrite('sound1.wav',x,8000);
sound(x,8000);