% Vincoli di missione:
%	Decollo da Sea Level con corsa 1800 m
%	Crociera a Ma = 0.85 a h = 10000 m
%	Virata in crociera a n = 1.1 (phi approx. 25°)
%	Autonomia 8000 km con 100 passeggeri

addpath('.\constraints_pkg');

global R gamma delta p_rif T_rif

R = 287;
gamma = 1.4;
delta = (gamma - 1) / 2;
[~, ~, T_rif, p_rif] = tropos(0);
g0 = 9.81;

AR = 11;
osw = 0.85;
k1 = 1 / (pi * AR * osw);
kto = 1.2;
CD0 = 0.02;
CL_max_to = 2.2;
CL_max_ln = 2.4;

carico_alare = linspace(3000, 9000, 1000);

% Vincolo: decollo
h = 0;
M0 = 0;
rho = tropos(h);

beta = 1;
alpha = alpha_break(0, 0);
corsa = 1800;

TWR_TO = beta ^ 2 / alpha * (kto ^ 2) / (corsa * rho * g0 * CL_max_to) * carico_alare;

% Vincolo: crociera
h = 10000;
M0 = 0.85;
[rho, a] = tropos(h);

beta = 0.9;
alpha = alpha_break(h, M0);
vel = M0 * a;
q = 0.5 * rho * vel ^ 2;
Ps = 0;
n = 1;

TWR_cr = beta / alpha * (q ./ (beta * carico_alare) .* (k1 .* (n * beta / q * carico_alare) .^ 2 + CD0) + Ps / vel);

% Vincolo: virata in crociera
h = 10000;
M0 = 0.85;
[rho, a] = tropos(h);

beta = 0.9;
alpha = alpha_break(h, M0);
vel = M0 * a;
q = 0.5 * rho * vel ^ 2;
Ps = 0;
n = 1.1;

TWR_vir = beta / alpha * (q ./ (beta * carico_alare) .* (k1 .* (n * beta / q * carico_alare) .^ 2 + CD0) + Ps / vel);

% Determinazione W_TO
% Da analisi di mercato
empty_over_to = 0.56;

% Per 100 passeggeri + bagaglio, da normativa
W_P = 9300 * g0;

% Per calcolo W_F (considerando N e non kg)
c1 = 0.0581e-4;
c2 = 0.25e-4;
M_cr = 0.85;
h_cr = 10000;
[~, a_cr] = tropos(h_cr);
range = 8000e3;

fuel_over_to = 1 - exp(-(c1 / M_cr + c2) * sqrt(4 * CD0 / (pi * AR * osw)) * range / a_cr * g0);

W_TO = W_P / (1 - fuel_over_to - empty_over_to);

err = 1;

W_TO_vec = [];

while err > 1e-2 
	W_TO_vec(end + 1) = W_TO;
	W_TO_old = W_TO;
	empty_over_to = 1.02 * W_TO ^ (-0.06);
	W_TO = W_P / (1 - fuel_over_to - empty_over_to);
	err = abs(W_TO - W_TO_old);
end
	
T_SL = W_TO * 0.3;

% Grafico Constraint Analysis
figure; hold on; grid on; leg = {};
plot(carico_alare, TWR_TO); leg{end + 1} = 'Decollo';
plot(carico_alare, TWR_cr); leg{end + 1} = 'Crociera';
plot(carico_alare, TWR_vir); leg{end + 1} = 'Virata';
plot(6000, 0.3, 'ko'); leg{end + 1} = 'Punto di progetto';
xlabel('W_{TO} / S'); xlim([3000, 9000]);
ylabel('T_{SL} / W_{TO}'); ylim([0, 0.6]);
legend(leg, 'location', 'northeast');

% saveTightFigure(gcf, 'constraint_analysis.pdf');

% Grafico convergenza W_TO
figure(2); hold on; grid on;
plot(W_TO_vec, '*');
xlabel('n_{step}');
ylabel('W_{TO} stimato [N]');
saveTightFigure(gcf, 'convergenza_peso.pdf');