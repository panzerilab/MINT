clc, clear,
% Definieren der Dimensionen
S = 10; % Anzahl der Proben
T = 5;  % Anzahl der Zeiteinheiten
N = 8;  % Anzahl der Neuronen
P = 3;  % Ein Parameter für die Dimension von B_tem und H
L = 4;  % Ein Parameter für die Dimension von B_spa
R = rand(N, T, S);
%[Acal,Wi,Wb,vaf,err] = sbtnmf(R,P,L,S);


[B_tem, H, B_spa, errorSum] = sNM3F(R, P, L);

R = reshape(R, [N, T*S]);
%[Acal,Wi,Wb,vaf,err] = sbtnmf(R,P,L,S);

disp(B_tem);