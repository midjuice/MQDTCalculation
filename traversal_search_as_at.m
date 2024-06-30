clear
clc

load("samples2.mat");
surf(X, Y, samples(:,:,end));
xlabel("Ks");
ylabel("Kt");
zlabel("Distance");