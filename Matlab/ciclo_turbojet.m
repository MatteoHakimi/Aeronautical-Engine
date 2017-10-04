addpath('.\ciclo_pkg');

% Parametri aria e carburante
global R gamma Cp
R = 287;
gamma = 1.4;
Cp = R * gamma / (gamma - 1);
Qf = 45 * 10^6;

T_req = 74.86e3;
mdot = 1;

fprintf('1) beta_c\n');
fprintf('2) T4 [K]\n');
fprintf('3) h [km]\n');
fprintf('4) Mach di crociera\n');
scelta = input('Inserisci la grandezza da variare > ');

i_min = input('Inserisci il valore minimo > ');
i_max = input('Inserisci il valore massimo > ');
di = (i_max - i_min) / 29;

j = 0;

h = 10000;
M = 0.8;
T4 = 1400;
beta_c = 12;

leg = {};
TSFC_vec = zeros(1, 30);
Ia_vec = zeros(1, 30);
eta_th_vec = zeros(1, 30);
eta_p_vec = zeros(1, 30);
eta_o_vec = zeros(1, 30);



for i = i_min : di : i_max
	j = j + 1;
	
	switch scelta
		case 1
			beta_c = i;
			xl = '\beta_c';
			xl_ad = xl;
		case 2
			T4 = i;
			xl = 'T_4 [K]';
			xl_ad = 'T_4';
		case 3
			h = i;
			xl = 'h [m]';
			xl_ad = 'h';
		case 4
			M = i;
			xl = 'M';
			xl_ad = 'M';
	end

	% Rendimenti
	eta_d = 0.95;
	eta_ac = 0.92;
	eta_mc = 0.96;
	eta_b = 0.9;
	eta_pb = 0.98;
	eta_at = 0.96;
	eta_mt = 0.99;
	eta_n = 0.94;

	% Condizioni 1: atmosfera standard
	[rho_a, a_a, T_a, p_a, ~, ~] = atmos(h);
	v0 = M * a_a;
	s_a = sF(p_a, T_a, R, gamma);

	% Presa dinamica
	[p2, T2, s2] = dinamica(p_a, T_a, eta_d, M);

	% Compressore
	[p3, T3, s3, Lc] = compressore(p2, T2, beta_c, eta_ac);

	% Combustore
	[p4, s4, f] = combustore(p3, T3, T4, eta_pb, eta_b, Qf);

	% Turbina
	[p5, T5, s5] = turbina(p4, T4, Lc, eta_at, eta_mt, eta_mc);

	% Ugello
	[T9, ue, s9] = ugello(p5, T5, p_a, eta_n);

	% Calcolo prestazioni
	[TSFC, Ia, eta_th, eta_p, eta_o, mdot] = prestazioni(ue, f, Qf, v0, T_req);
	TSFC_vec(j) = TSFC;
	Ia_vec(j) = Ia;
	eta_th_vec(j) = eta_th;
	eta_p_vec(j) = eta_p;
	eta_o_vec(j) = eta_o;

	% Grafico ciclo
	plot_ciclo;
% 	legend(leg, 'location', 'northwest');
end
% saveTightFigure(gcf, 'ciclo.pdf');

figure(2);
plot(i_min : di : i_max, TSFC_vec);
xlabel(xl);
ylabel('TSFC [kg / N h]');
title('Variazione TSFC');
% saveTightFigure(gcf, 'TSFC.pdf');

figure(3)
plot(i_min : di : i_max, Ia_vec);
xlabel(xl);
ylabel('I_a [m/s]');
title('Variazione I_a');
% saveTightFigure(gcf, 'Ia.pdf');

figure(4)
plot(i_min : di : i_max, eta_th_vec);
xlabel(xl);
ylabel('\eta_{th}');
title('Variazione \eta_{th}');

figure(5)
plot(i_min : di : i_max, eta_p_vec);
xlabel(xl);
ylabel('\eta_p');
title('Variazione \eta_p');

figure(6)
plot(i_min : di : i_max, eta_o_vec);
xlabel(xl);
ylabel('\eta_o');
title('Variazione \eta_o');

if scelta == 4
	coef = polyfit(i_min : di : i_max, TSFC_vec, 1);
end