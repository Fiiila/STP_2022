close  all
clear all
clc

%% Priklad 1
%% TEORETICKY VYPOCET
TAU = 0:1:5;
K = 0:1:94;
N = 100; %pocet casovych okamziku
M = 10e4; % pocet realizaci procesu


%BILY SUM Wk
%VAR_W = ones(N,1); % variance bileho sumu
%E_W = ones(N,1);   % stredni hodnoty bileho sumu
Q = 3;
b = 0.5;
T = 1;

%stredni hodnota
E_W = 0;
%E_W = E_W*E;
%variance (rozptyl)
VAR_W = Q*(1-exp(-2*b*T));
%VAR_W = VAR_W*VAR;

%proces Xk
VAR_X0 = Q;
E_X0 = 0;
E_X = zeros(N,1);
E_X(1) = E_X0; % dano pocatecni podminkou
for k=2:1:N
    E_X(k) = exp(-b*T)*E_X(k-1)+E_W;
end

% teoreticke kovariancni funkce
teor_kov_fce = zeros(length(TAU),1);
for tau=1:1:length(TAU)
    teor_kov_fce(tau) = Q*exp(-TAU(tau)*b*T); % zjistena funkce manualnim vypoctem
end

%% GENEROVANI REALIZACI PROCESU

X = zeros(N,M);

for i=1:1:M
    X(1,i) = randn* sqrt(VAR_X0);% dano pocatecni podminkou
    for k=1:1:N
        X(k+1,i) = exp(-b*T)*X(k,i) + randn * sqrt(VAR_W);
    end
end

figure
plot(0:1:N,X(:,1:10))
xlim([0,100])
title('10 realizací zadaného procesu')
xlabel('krok [k]')
ylabel('hodnota [X]')

set(gcf,'Renderer','Painter');
saveas(gcf,'priklad_01_realizace.png','png')
saveas(gcf,'priklad_01_realizace.eps','epsc')
%% ODHAD KOVARIANCE

odhad_kov_fce = zeros(length(K),length(TAU));
for k=1:1:length(K)
    for tau=1:1:length(TAU)
        Xk = X(k,:);
        Xk_tau = X(k+tau-1,:);
        E_Xk = mean(Xk);
        E_Xk_tau = mean(Xk_tau);
        for n = 1:1:M
            odhad_kov_fce(k,tau) = odhad_kov_fce(k,tau) + (Xk(n)-E_Xk)*(Xk_tau(n)-E_Xk_tau);
        end
        odhad_kov_fce(k,tau) = odhad_kov_fce(k,tau)/M;
    end
end

plot(0:1:K(end),odhad_kov_fce, 'LineWidth', 2)
hold on
plot(0:1:K(end),ones(length(K),length(TAU)).*teor_kov_fce', 'black--')
legend('COV[X_{k},X_{k}]',...
       'COV[X_{k},X_{k+1}]',...
       'COV[X_{k},X_{k+2}]',...
       'COV[X_{k},X_{k+3}]',...
       'COV[X_{k},X_{k+4}]',...
       'COV[X_{k},X_{k+5}]',...
       'teoreticky vypočtené \newlinehodnoty COV[X_{k},X_{k+\tau}]')
title("Porovnání vypočtených a odhadnutých kovariancí")
xlabel("krok [k]")
ylabel("hodnota COV[X_{k},X_{k+\tau}]")
saveas(gcf,'priklad_01_porovnani_kovarianci.png','png')
saveas(gcf,'priklad_01_porovnani_kovarianci.eps','epsc')
%% Priklad 2
%% generovani realizaci Wienerova procesu

X0 = 0;
E_W = 0;
VAR_W = 1;
TAU = 0:1:5;
K = 0:1:94;
M = 10e4;
N = 100;

X = zeros(N,M);
for i=1:1:M
    X(1,i)=X0;
    for k=1:1:N
        X(k+1,i) = X(k,i) + randn*sqrt(VAR_W);
    end
end
figure
plot(0:1:N,X(:,1:8))
xlim([0,100])
title('8 realizací zadaného Wienerova procesu')
xlabel('krok [k]')
ylabel('hodnota [X]')

saveas(gcf,'priklad_02_realizace_procesu.png','png')
saveas(gcf,'priklad_02_realizace_procesu.eps','epsc')

%% teoreticky vypocet kovariancni funkce

%% odhad kovariancni funkce
odhad_kov_fce = zeros(length(K),length(TAU));
for k=1:1:length(K)
    for tau=1:1:length(TAU)
        Xk = X(k,:);
        Xk_tau = X(k+tau-1,:);
        E_Xk = mean(Xk);
        E_Xk_tau = mean(Xk_tau);
        for n = 1:1:M
            odhad_kov_fce(k,tau) = odhad_kov_fce(k,tau) + (Xk(n)-E_Xk)*(Xk_tau(n)-E_Xk_tau);
        end
        odhad_kov_fce(k,tau) = odhad_kov_fce(k,tau)/M;
    end
end

plot(K,odhad_kov_fce, 'LineWidth', 2)
hold on
grid on
plot(K,K,'black--')
legend('COV[X_{k},X_{k}]',...
       'COV[X_{k},X_{k+1}]',...
       'COV[X_{k},X_{k+2}]',...
       'COV[X_{k},X_{k+3}]',...
       'COV[X_{k},X_{k+4}]',...
       'COV[X_{k},X_{k+5}]',...
       'teoreticky vypočtené \newlinehodnoty COV[X_{k},X_{k+\tau}]',...
       'Location','SouthEast')
title("Porovnání vypočtených a odhadnutých kovariancí")
xlabel("krok [k]")
ylabel("hodnota COV[X_{k},X_{k+\tau}]")
saveas(gcf,'priklad_02_porovnani_kovarianci.png','png')
saveas(gcf,'priklad_02_porovnani_kovarianci.eps','epsc')

%% Priklad 3
%% generovani realizaci systemu
E_W = 0;
VAR_W = 3;
E_V = 0;
VAR_V = 2;
E_X0 = 1;
VAR_X0 = 5;
K = 0:1:99;

X = zeros(N,M);
Z = zeros(N,M);
for i=1:1:M
    X(1,i) = E_X0 + randn*sqrt(VAR_X0);
    Z(1,i) = 5*X(1,i) + randn*sqrt(VAR_V);
    for k=1:1:N-1
        X(k+1,i) = 0.95*X(k,i) + 0.5*randn*sqrt(VAR_W);
        Z(k+1,i) = 5*X(k+1,i) + randn*sqrt(VAR_V);
    end
end
% vykresleni X
figure
plot(0:1:N-1,X(:,1:10))
xlim([0,100])
title('10 realizací X_k zadaného systému')
xlabel('krok [k]')
ylabel('hodnota [X]')

saveas(gcf,'priklad_03_realizace_procesu_X.png','png')
saveas(gcf,'priklad_03_realizace_procesu_X.eps','epsc')

%vykresleni Z
figure
plot(0:1:N-1,Z(:,1:10))
xlim([0,100])
title('10 realizací Z_k zadaného systému')
xlabel('krok [k]')
ylabel('hodnota [Z]')

saveas(gcf,'priklad_03_realizace_procesu_Z.png','png')
saveas(gcf,'priklad_03_realizace_procesu_Z.eps','epsc')

%%odhady

odhad_E_X = zeros(length(K),1);
odhad_E_Z = zeros(length(K),1);
odhad_VAR_X = zeros(length(K),1);
odhad_VAR_Z = zeros(length(K),1);
for k=1:1:length(K)
    odhad_E_X(k) = mean(X(k,:));
    odhad_E_Z(k) = mean(Z(k,:));
    for n = 1:1:M
        odhad_VAR_X(k) = odhad_VAR_X(k)+(X(k,n)-odhad_E_X(k))^2;
        odhad_VAR_Z(k) = odhad_VAR_Z(k)+(Z(k,n)-odhad_E_Z(k))^2;
    end
    odhad_VAR_X(k) = odhad_VAR_X(k)/M;
    odhad_VAR_Z(k) = odhad_VAR_Z(k)/M;
end
%vykresleni strednich hodnot X
figure
x = 0:1:N-1;
plot(x,odhad_E_X)
hold on
y = 0.95.^(x)
plot(x,y,'black--')
legend('odhad E[X_k]','vypočtené E[X_k]')
title('Porovnání odhadnutých a vypočtených hodnot E[X_k]')
xlabel("krok [k]")
ylabel("hodnota E[X_k]")
saveas(gcf,'priklad_03_porovnani_E_X.png','png')
saveas(gcf,'priklad_03_porovnani_E_X.eps','epsc')

%vykresleni strednich hodnot Z
figure
x = 0:1:N-1;
plot(x,odhad_E_Z)
hold on
y = 5*0.95.^(x)
plot(x,y,'black--')
legend('odhad E[Z_k]','vypočtené E[Z_k]')
title('Porovnání odhadnutých a vypočtených hodnot E[Z_k]')
xlabel("krok [k]")
ylabel("hodnota E[Z_k]")
saveas(gcf,'priklad_03_porovnani_E_Z.png','png')
saveas(gcf,'priklad_03_porovnani_E_Z.eps','epsc')

%vykresleni varianci X
figure
x = 0:1:N-1;
plot(x,odhad_VAR_X)
hold on
sumas=zeros(1,length(x));
for i=1:1:length(x)
    for j=1:1:i-1
        sumas(i) = sumas(i) + 0.95^(2*x(j));
    end
end
y = 0.95.^(2*x)*5+0.5^2*3.*sumas;
teor_VAR_X = y;
plot(x,y,'black--')
legend('odhad VAR[X_k]','vypočtené VAR[X_k]','Location','SouthEast')
title('Porovnání odhadnutých a vypočtených hodnot VAR[X_k]')
xlabel("krok [k]")
ylabel("hodnota VAR[X_k]")
saveas(gcf,'priklad_03_porovnani_VAR_X.png','png')
saveas(gcf,'priklad_03_porovnani_VAR_X.eps','epsc')

%vykresleni varianci Z
figure
x = 0:1:N-1;
plot(x,odhad_VAR_Z)
hold on
y = 5^2*teor_VAR_X+VAR_V;
plot(x,y,'black--')
legend('odhad VAR[Z_k]','vypočtené VAR[Z_k]','Location','SouthEast')
title('Porovnání odhadnutých a vypočtených hodnot VAR[Z_k]')
xlabel("krok [k]")
ylabel("hodnota VAR[Z_k]")
saveas(gcf,'priklad_03_porovnani_VAR_Z.png','png')
saveas(gcf,'priklad_03_porovnani_VAR_Z.eps','epsc')

