N = 500;
x = linspace(-20,20,N);
y = zeros(1, N);
for index=1:N
    y(index) = Lambdal(x(index),-100,2,10);
end

plot(x,y,'o')
ylim([-100,100])