clc; clear; close all;

% Szimbólumok deklarálása
syms p0 p1 p2 p us uc il duc dil %p = p- = p+
syms L C R

% Egyenletrendszer felírása
eq1 = (C*duc) + (uc-p1)/R == 0;
eq2 = (p1-uc)/R +(p1-p0)/(2*R) -il == 0;
eq3 = il+ (p1+(L*dil)-p2)/(2*R) == 0;
eq4 = (p2-p1+(L*dil))/(2*R)+(p2-p)/R +(p2-us)/R== 0;
eq5 = p/R+(p-p0)/R == 0;
eq6 = (p-p2)/R+(p-p0)/(0.1*R) == 0;

% Az egyenletek listája
eqs = [eq1, eq2, eq3, eq4, eq5, eq6];
%kell : p2 dil duc
% Megoldandó változók
vars = [p0,p1,p2,p,uc,il];

% Egyenletrendszer megoldása
sol = solve(eqs, vars);

% derivált kifejezések visszaállítása (duc, dil)
feq1 = uc == sol.uc;%rendezett e.gy.
feq2 = il == sol.il;%rendezett e.gy.
fsol = solve([feq1, feq2], [duc, dil]);

% Helyettesítjük a deriváltakat a p1 kifejezésébe
p1_expr = subs(sol.p1, [duc, dil], [fsol.duc, fsol.dil]);
dil_expr = fsol.dil;

% Kiszámítjuk a p1 + L*dil kifejezést
p1_plus_Ldil = p1_expr + L*dil_expr;


% Helyettesítjük a deriváltakat az u1 kifejezésébe
p2_expr = subs(sol.p2, [duc, dil], [fsol.duc, fsol.dil]);
% Kiíratás
disp("duc = ");
disp(collect(expand(fsol.duc), [uc, il, us]));

disp("dil = ");
disp(collect(expand(fsol.dil), [uc, il, us]));


disp("p1 + L*dil = ");
disp(collect(expand(p1_plus_Ldil), [uc, il, us]));


%behelyettesitett rendszer matrixok
A = [-15/(43*21*4),27/(43*21);-19/(43*2200),-(181*4)/(86*2200)];
B = [-1/(43*4*21);29/(86*2200)];
C = [9/43;-127*4/86];
D = 27/86;%us
disp("A matrix:");
disp(A);

%% 2.2
disp("sajatertekek:");
sajat_ert = eig(A);
lambda_1 = sajat_ert(1);
lambda_2 = sajat_ert(2);
disp(lambda_1);
disp(lambda_2);

%% 2.3
m1=[1;(real(lambda_1)-A(1,1))/A(1,2)];%sajatvektor 1
m2=[1;(real(lambda_2)-A(2,1))/A(2,2)];%sajatvektor 2
M = [m1,m2];% "sajat matrix"

%gerjesztes valsz 
xg = -inv(A)*B;

%k1 k2 szamitasa 
k = -inv(M)*xg;
k1 = k(1);
k2 = k(2);

g1=(C'*k1*m1);
g2=(C'*k2*m2);
g3 = ((C'*xg)+D);

% Időtartomány definiálása (optimalizált a rendszer dinamikája alapján)
time_constant = max(abs(1./real(sajat_ert)));
t = linspace(0, 5*time_constant, 1000); % 5 időállandóig számolunk

% Ugrásválasz számítása
g_t = g1*exp(lambda_1*t) + g2*exp(lambda_2*t) + g3;

% Impulzusválasz (súlyfüggvény) - az ugrásválasz deriváltja
h_t = g1*lambda_1*exp(lambda_1*t) + g2*lambda_2*exp(lambda_2*t);

%% Ábrázolás
figure('Position', [100, 100, 900, 600]);

% 1. Ugrásválasz
subplot(2,1,1);
plot(t, g_t, 'b', 'LineWidth', 2);
hold on;
plot([t(1) t(end)], [g3 g3], 'r--', 'LineWidth', 1.5);
grid on;
xlabel('Idő (t) [s]', 'FontSize', 10);
ylabel('Amplitúdó', 'FontSize', 10);
title('Ugrásválasz (átmeneti függvény)', 'FontSize', 12);
legend('g(t)', 'Stacionárius állapot', 'Location', 'best');
xlim([t(1) t(end)]);

% 2. Impulzusválasz
subplot(2,1,2);
plot(t, h_t, 'Color', [0 0.5 0], 'LineWidth', 2);
grid on;
xlabel('Idő (t) [s]', 'FontSize', 10);
ylabel('Amplitúdó', 'FontSize', 10);
title('Impulzusválasz (súlyfüggvény)', 'FontSize', 12);
legend('h(t)', 'Location', 'best');
xlim([t(1) t(end)]);

%% Stabilitás elemzés
figure;
plot(real(sajat_ert), imag(sajat_ert), 'x', 'MarkerSize', 12, 'LineWidth', 2);
hold on;
plot(0, 0, 'r+', 'LineWidth', 2);
grid on;
xlabel('Re(\lambda)', 'FontSize', 10);
ylabel('Im(\lambda)', 'FontSize', 10);
title('Sajátértékek a komplex síkon', 'FontSize', 12);
legend('Sajátértékek', 'Origó', 'Location', 'best');

% Stabilitás megállapítása
if all(real(sajat_ert) < 0)
    stability = 'STABILIS';
    color = [0 0.7 0];
else
    stability = 'INSTABILIS';
    color = [0.9 0 0];
end

annotation('textbox', [0.2, 0.05, 0.6, 0.05], 'String',...
    ['Rendszer állapota: ' stability ' (λ1 = ' num2str(lambda_1) ', λ2 = ' num2str(lambda_2) ')'],...
    'FitBoxToText', 'on', 'BackgroundColor', color,...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 10);

%% Konzol eredmények
disp('=== Eredmények ===');
disp(['Sajátérték 1 (λ1): ' num2str(lambda_1)]);
disp(['Sajátérték 2 (λ2): ' num2str(lambda_2)]);
disp(['Stacionárius érték (g3): ' num2str(g3)]);
disp(['Rendszer: ' stability]);
disp('=================');
