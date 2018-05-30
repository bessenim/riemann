x = linspace(-0.5, 0.5, 1001);
u1 = zeros(1001,1);
u2 = zeros(1001,1);
u3 = zeros(1001,1);
u4 = zeros(1001,1);
u5 = zeros(1001,1);
for i=1:1001
    if i <= 100
        u1(i) = -4;
        u2(i) = -3;
        u3(i) = 1;
        u4(i) = -2;
        u5(i) = 1;
    end
    if (100<i && i<=300)
        u1(i) = 0;
        u2(i) = -1;
        u3(i) = 1;
        u4(i) = 0;
        u5(i) = 1;
    end
    if (300<i && i<=500)
        u1(i) = 0;
        u2(i) = -1;
        u3(i) = 0;
        u4(i) = 0;
        u5(i) = 0;
    end
    if (500<i && i<=700)
        u1(i) = 0;
        u2(i) = 0;
        u3(i) = 0;
        u4(i) = 0;
        u5(i) = 0;
    end
    if (700<i && i<=900)
        u1(i) = 0;
        u2(i) = 0;
        u3(i) = 2;
        u4(i) = 0;
        u5(i) = -2;
    end
    if i>900
        u1(i) = -4;
        u2(i) = -2;
        u3(i) = 2;
        u4(i) = 2;
        u5(i) = -2;
    end
end
plot(x, u1, 'DisplayName','sigma11')
hold on
plot(x, u2, 'DisplayName','sigma22')
hold on
plot(x, u3, 'DisplayName','sigma12')
hold on
plot(x, u4, 'DisplayName','v1')
hold on
plot(x, u5, 'DisplayName','v2')
hold off
ylim([-5 3])

legend({'sigma11', 'sigma22', 'sigma12', 'v1', 'v2'}, 'Location', 'north')

    