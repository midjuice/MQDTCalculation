clc
clear


Kc = zeros(5,5,3);
de = zeros(1,3);
chik = zeros(4,4,3);
ll = [1,3,5];
for i = 1:length(ll)
    [Kc(:,:,i), ~, chik(:,:,i), ~, ~, ~, ~,de(i)] = Kcfunc(0, 'Li6','ab',100,ll(i),2.003);
end

chis = zeros(4,3);

for i =1:3
    chis(:,i)= squeeze(diag(chik(:,:,i)));
end

kc11 = squeeze(Kc(1,1,:));
kc13 = squeeze(Kc(1,3,:));
kc24 = squeeze(Kc(2,4,:));
kc45 = squeeze(Kc(4,5,:));


    
plot(ll,kc11,'-o','LineWidth', 2);
hold on;
plot(ll,kc13,'-o','LineWidth', 2);
hold on;
plot(ll,kc24,'-o','LineWidth', 2);
hold on;
plot(ll,kc45,'-o','LineWidth', 2);


legend("K_{11}", "K_{13}", "K_{24}" ,"K_{45}");
ylabel("K");
xlabel("Partical wave $\mathcal{L}$",Interpreter="latex");
set(gca, 'XTick', [1,3,5]);