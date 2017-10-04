addpath('.\compressore_pkg');

clear all 
close all 
clc

global R gamma C_P delta
R = 287;
gamma = 1.4;
C_P = R * gamma / (gamma - 1);
delta = 0.2;

p_in = 28509.33;
T_in = 231.652;
M_in = 0.45;
M_tip_in = 0.9;
beta_tot = 15;
chi = 0.7;
AR = 1.8;
DF = 0.45;
m_dot = 36.15;
eta_p = 0.88;
lam_rm_in = 0.5;

[p1_vec, T1_vec, p2_vec, T2_vec, p3_vec, T3_vec, M2, M3, lam_vec, num_stadi, num_pale, xi_rot, xi_stat, h_vec, r_tip_in] = ...
	comp_sizing(p_in, T_in, M_in, M_tip_in, beta_tot, chi, AR, DF, m_dot, eta_p, lam_rm_in);

r_hub = r_tip_in - h_vec ; 
% Pressione
figure(1);
press_matrix = [p2_vec; p3_vec];
press_vec = [];
press_lab = {'0'};
for i = 1 : num_stadi
	press_vec(end + 1) = press_matrix(1, i);
	press_vec(end + 1) = press_matrix(2, i);
	press_lab{end + 1} = sprintf('%dr', i);
	press_lab{end + 1} = sprintf('%ds', i);
end
n_stadi=1:num_stadi ;
r_tip_vec = r_tip_in*ones (1 , num_stadi); 

press_vec = [p1_vec(1), press_vec];
bar(press_vec);
set(gca, 'xtick', 1 : 2 * num_stadi + 1, 'xticklabel', press_lab);
xlabel('Numero stadio');
ylabel('Pressione [Pa]');
title('Andamento della pressione in uscita dalle schiere')
%saveTightFigure(gcf, 'pressione_uscita.pdf');

% Temperatura
figure(2);
temp_matrix = [T2_vec; T3_vec];
temp_vec = [];
for i = 1 : num_stadi
	temp_vec(end + 1) = temp_matrix(1, i);
	temp_vec(end + 1) = temp_matrix(2, i);
end

temp_vec = [T1_vec(1), temp_vec];
bar(temp_vec);
set(gca, 'xtick', 1 : 2 * num_stadi + 1, 'xticklabel', press_lab);
xlabel('Numero stadio');
ylabel('Temperatura [K]');
title('Andamento della temperatura in uscita dalle schiere')
%saveTightFigure(gcf, 'temperatura_uscita.pdf');

% Mach
figure(4);
mach_matrix = [M2; M3];
mach_vec = [];
for i = 1 : num_stadi
	mach_vec(end + 1) = mach_matrix(1, i);
	mach_vec(end + 1) = mach_matrix(2, i);
end

mach_vec = [M_in, mach_vec];
bar(mach_vec);
set(gca, 'xtick', 1 : 2 * num_stadi + 1, 'xticklabel', press_lab);
bar(M3);
xlabel('Numero stadio');
ylabel('Ma');
title('Andamento del Mach in uscita dagli stadi')
%saveTightFigure(gcf, 'mach_uscita.pdf');

% Grado di reazione
figure(3);
bar(lam_vec);
xlabel('Numero stadio');
ylabel('\Lambda_m');
title('Andamento grado di reazione negli stadi, raggio medio');
%saveTightFigure(gcf, 'lambda_stadi.pdf');

% Perdite rotore
figure(5);
bar(xi_rot);
xlabel('Numero di stadio');
ylabel('\xi_{rot}');
title('Andamento delle perdite di entalpia nei rotori');
%saveTightFigure(gcf, 'xi_rot.pdf');

% Perdite statore
figure(6);
bar(xi_stat);
xlabel('Numero di stadio');
ylabel('\xi_{stat}');
title('Andamento delle perdite di entalpia negli statori');
%saveTightFigure(gcf, 'xi_stat.pdf');

% Distribuzione raggio hub
figure(7);
plot ( n_stadi, r_hub,'b')
hold on
plot ( n_stadi, r_tip_vec,'b' )
xlabel('Numero di stadio');
ylabel('r_{hub} [ m ]');
axis( [ 1 num_stadi 0 0.7 ])
title('Distribuzione raggio hub ');
%saveTightFigure(gcf, 'xi_stat.pdf');