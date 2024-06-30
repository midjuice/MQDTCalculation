clear
clc

rsize = 1000;
rs = linspace(0.01,100,rsize);
Bs = linspace(.1,1100,10);
prop = zeros(5,1);
for i = 1:length(Bs)
    [~,ratio] = norEfun(0,'Li6','ab',Bs(i),2,2.003);
    prop(:,i) = ratio;

end


for i = 1:5
    plot(Bs,prop(i,:),'-o', 'MarkerSize', 4,'LineWidth', 1.5);
    hold on;
end
legend("aa", "bh", "fh" ,"ag", "gg");
ylabel("Proportion");
xlabel("B(G)");
% rect_position = [1030, 0.14, 30, 0.12];
% rectangle('Position', rect_position, 'EdgeColor', 'r', 'LineWidth', 1, 'FaceColor', 'none', 'LineStyle', '--');
% disp(1);
% 
% axes('Position', [0.55, 0.15, 0.35, 0.29]);
% box on; % 使嵌入图有边框
% prop1 = zeros(5,1);
% % 在嵌入坐标轴中绘制细节图
% Bs1 = [1038,1039,1045,1050,1054,1055];
% for i = 1:length(Bs1)
%     [~,ratio] = norEfun(0,'Li7','aa',Bs1(i),2,2.003);
%     prop1(:,i) = ratio;
% 
% end
% 
% for i = 1:5
%     plot(Bs1,prop1(i,:),'-o', 'MarkerSize', 2,'LineWidth', 0.75);
%     hold on;
% end
% 
% ylabel("Proportion");
% xlabel("B(G)");


