clc
clear
size = 1000;
Bs = linspace(.1,1500,size);
Kc = zeros(5,5,size);
de = zeros(1,size);
chik = zeros(4,4,size);

for i = 1:size
    [Kc(:,:,i), ~, chik(:,:,i), ~, ~, ~, ~,de(i)] = Kcfunc(0, 'Li7','aa',Bs(i),2,2.003);
end
chis = zeros(4,size);

for i =1:size
    chis(:,i)= squeeze(diag(chik(:,:,i)));
end
% 
% kc11 = squeeze(Kc(1,1,:));
% kc13 = squeeze(Kc(1,3,:));
% kc24 = squeeze(Kc(2,4,:));
% kc45 = squeeze(Kc(4,5,:));


    
% plot(Bs,kc11,'LineWidth', 2);
% hold on;
% plot(Bs,kc13,'LineWidth', 2);
% hold on;
% plot(Bs,kc24,'LineWidth', 2);
% hold on;
% plot(Bs,kc45,'LineWidth', 2);


% legend("K_{11}", "K_{13}", "K_{24}" ,"K_{45}");
% ylabel("K");
% xlabel("B(G)");

% 
% for i = 1:4
%     plot(Bs,chis(i,:),'LineWidth', 2);
%     hold on;
% end
% legend('\chi_{11}','\chi_{22}','\chi_{33}','\chi_{44}');
% xlabel('B(G)');
% ylabel('\chi_{\Delta=0}');

plot(Bs,de,'LineWidth', 2);
xlabel('B(G)');
ylabel('|\chi-K_c|');

rect_position = [1035, -0.1, 25, 0.2];
rectangle('Position', rect_position, 'EdgeColor', 'r', 'LineWidth', 1, 'FaceColor', 'none', 'LineStyle', '--');



axes('Position', [0.55, 0.15, 0.35, 0.29]);
box on; % 使嵌入图有边框
prop1 = zeros(5,1);
% 在嵌入坐标轴中绘制细节图
Bs1 = linspace(1035.1,1060,200);
size1 = length(Bs1);
de1 = zeros(1,size1);
for i = 1:size1
    [~, ~, ~, ~, ~, ~, ~,de1(i)] = Kcfunc(0, 'Li7','aa',Bs1(i),2,2.003);
end

plot(Bs1,de1,'-','LineWidth', 0.75);
yline(0,'--r');
xlabel('B(G)');
ylabel('|\chi-K_c|');




