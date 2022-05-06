z = 0:0.21:400;
t = sqrt(z);
t = [-t, t];
t = sort(t);
x = fresnelc(t);
y = fresnels(t);
figure; hold on; grid on; plot(x, y);
