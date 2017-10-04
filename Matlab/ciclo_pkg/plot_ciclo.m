figure(1); hold on;
[s_d, T_d] = ramo_entro(p2, p_a, T2, T_a, s2, s_a);
[s_c, T_c] = ramo_entro(p3, p2, T3, T2, s3, s2);
[s_b, T_b] = ramo_entro(p4, p3, T4, T3, s4, s3);
[s_t, T_t] = ramo_entro(p5, p4, T5, T4, s5, s4);  
[s_n, T_n] = ramo_entro(p_a, p5, T9, T5, s9, s5);
[s_r, T_r] = ramo_entro(p_a, p_a, T_a, T9, s_a, s9);

% colors = get(gca, 'colororder');
% colors = [colors; colors; colors; colors; colors];

r = [1, 0, 0];
b = [0, 0, 1];

colors = b + (r - b) * (j - 1) / 30;

clear s_plot T_plot

s_plot = [s_d, s_c, s_b, s_t, s_n, s_r];
T_plot = [T_d, T_c, T_b, T_t, T_n, T_r];

linee = plot(s_plot, T_plot);
set(linee, 'color', colors);
leg{end + 1} = sprintf('%s = %6.2f', xl_ad, i);

xlabel('s [J / kg K]');
ylabel('T [K]');
title('Ciclo Turbojet');