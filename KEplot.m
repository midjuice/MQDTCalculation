clc
clear

Es = linspace(-250,250,500);
size = length(Es);
Kc = zeros(5,5,size);
de = zeros(1,size);
chik = zeros(4,4,size);

for i = 1:size
    [Kc(:,:,i), ~, chik(:,:,i), ~, ~, ~, ~,de(i)] = Kcfunc(Es(i), 'Li6','ab',100,1,2.003);
end

chis = zeros(4,size);

for i =1:size
    chis(:,i)= squeeze(diag(chik(:,:,i)));
end

kc11 = squeeze(Kc(1,1,:));
kc13 = squeeze(Kc(1,3,:));
kc24 = squeeze(Kc(2,4,:));
kc45 = squeeze(Kc(4,5,:));


    
plot(Es,kc11,'LineWidth', 2);
hold on;
plot(Es,kc13,'LineWidth', 2);
hold on;
plot(Es,kc24,'LineWidth', 2);
hold on;
plot(Es,kc45,'LineWidth', 2);


legend("K_{11}", "K_{13}", "K_{24}" ,"K_{45}");
ylabel("K");
xlabel("$\varepsilon$",Interpreter="latex");


% for i = 1:4
%     plot(Bs,chis(i,:),'LineWidth', 2);
%     hold on;
% end
% legend('\chi_{11}','\chi_{22}','\chi_{33}','\chi_{44}');
% xlabel('B(G)');
% ylabel('\chi_{\Delta=0}');