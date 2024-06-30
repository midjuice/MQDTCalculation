clear
clc

points = 1000;
Es_one_over_three = linspace(-6, -0.1, points);
chi_array = zeros(1,1);

for i=1:length(Es_one_over_three)
    chi_array(i) = real(chilfunc(Es_one_over_three(i)^3, 1));
end

plot(Es_one_over_three, chi_array, '-o', markersize=1.5);
xlabel('$\varepsilon^{1/3}$','interpreter','latex');
ylabel('$\chi^c$','interpreter','latex');
