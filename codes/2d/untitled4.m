delta = 2.067142159891546e-04;
delta2 = delta/4.7;
delta3 = delta2/4.5;
delta4 = delta3/4.1;
delta5 = delta4/4;
x = [2,4,8,16,32];
x1 = log(x)
y = [delta, delta2, delta3, delta4, delta5];
y1 = log(y)
plot(x1,y1, ':r*')
xlabel('Логарифм степени дробления сетки')
ylabel('Логарифм максимального квадрата разности решений')