clc
clearvars
[signal, frekvence_vzorkovani] = audioread('/Users/patrikbuben/Desktop/BPC-SAS/270566.wav');
audiowrite('/Users/patrikbuben/Desktop/BPC-SAS/270566.wav',signal, frekvence_vzorkovani);
%sound(signal, frekvence_vzorkovani)

%zakladni parametry jako delka a frekvence
N = length(signal); %celkovy pocet vzorku
t = (0:N-1)/frekvence_vzorkovani; %casova osa od 0 do delky signalu [vzorec T = 1/f proste]

%transformace
F = fft(signal) / N; %fourier fast transformation
f = (0:N-1) * (frekvence_vzorkovani / N);

F_abs = abs(F);
f_plot = f(1:floor(N/2));
F_abs_plot = F_abs(1:floor(N/2));

figure('Name', 'Úkol 1: Filtrace signalu');
subplot(2, 1, 1)
plot(t, signal)
title('Puvodni signal (ptacci)')
xlabel('Cas [s]')
ylabel('Amplituda')
grid
hold on

subplot(2, 1, 2)
loglog(f_plot(2: end), F_abs_plot(2:end));
title('Puvodni amplitudove spektrum (obraz signalu)');
xlabel('Frekvence [Hz]')
ylabel('|F(f)|')
grid
hold on

F_filtrovana = F; %zkopiroval jsem si spektrum, abych ho mohl filtrovat
%ted podle grafu najdeme rusive frekvence ja jsem nasel (peaks kde to
%vystreli jak blazen)

f_peak = 10884; %[Hz]
sirka_pasma = 80; %timto odrezavame

for peak_freq = f_peak

    %tip4 - symetricka frekvence = spektrum real sig musi byt sym (je
    %zrcadlem --> "gumujeme" jak nas peak tak zrcadlovy
    sym_freq = frekvence_vzorkovani - peak_freq;

   %indexy pro hlavní frekvenci
    indexy_k_hlavni = find(f > (peak_freq - sirka_pasma/2) & f < (peak_freq + sirka_pasma/2));
                       
    %indexy pro symetrickou frekvenci, toto je pry ze zakona a z tipu4
    indexy_k_sym = find(f > (sym_freq - sirka_pasma/2) & f < (sym_freq + sirka_pasma/2));

    %vynulujeme hodnoty ve spektru
    F_filtrovana(indexy_k_hlavni) = 0;
    F_filtrovana(indexy_k_sym) = 0;
    
    fprintf('Filtruji pásmo kolem %.1f Hz a %.1f Hz\n', peak_freq, sym_freq);
end

f_filtrovana = ifft(F_filtrovana);
%chceme jenom real, tip4
f_filtrovana = real(f_filtrovana);

subplot(2, 1, 1)
plot(t, f_filtrovana, 'r');
legend('Puvodni signal', 'FIltrovany signal');
hold off;

F_filtrovana_abs = abs(F_filtrovana);
F_filtrovana_abs_plot = F_filtrovana_abs(1:floor(N/2));

subplot(2, 1, 2);
loglog(f_plot(2:end), F_filtrovana_abs_plot(2:end), "yellow");
legend('Puvodni spektrum', 'Filtrovane spektrum');
hold off;

%Tip 3
max_val = max(abs(f_filtrovana));
f_normal = f_filtrovana / max_val;

audiowrite('filtrovano.wav', f_normal, frekvence_vzorkovani);
[novy_signal, nova_fs] = audioread('/Users/patrikbuben/Desktop/BPC-SAS/cvikahome/filtrovano.wav');
sound(novy_signal, nova_fs);

%% filtr 2 radu (ukol 2)
p = tf('p');
Fs = frekvence_vzorkovani; 

f_ruseni = 10884;
w_ruseni = 2 * pi * f_ruseni; 

w0 = w_ruseni/3; 

%F(p) pro DPF: F(p) = w0^2 / (p^2 + 2*w0*p + w0^2)
system_P = (w0^2) / (p^2 + 2 * w0 * p + w0^2);

figure('Name', 'Úkol 2: Bodeho diagram Dorní propusti');
bode(system_P);
grid on;
title('Frekvenční charakteristika DPF');

%% diskretizace
z = tf('z', 1/frekvence_vzorkovani); 
Ts = 1 / frekvence_vzorkovani;

E = exp(-w0 * Ts);

% Koeficienty odvozené metodou ZOH:
num_b1 = 1 - E - w0*Ts*E;
num_b0 = E^2 - E + w0*Ts*E;
den_a1 = -2*E;
den_a0 = E^2;

system_Z = (num_b1 * z + num_b0) / (z^2 + den_a1 * z + den_a0);
system_Z = minreal(system_Z, 1e-4);


system_P_c2d = c2d(system_P, Ts, 'zoh'); 
%kdyz pouziji 'zoh' tak se kryji

figure('Name', 'Úkol 3: Diskretizace');

%Přechodová charakteristika
subplot(2, 2, 1);
step(system_P, 'b', system_Z, 'r--', system_P_c2d, 'g:'); 
legend('Spojitý', 'Analytický', 'Matlab c2d');
title('Přechodová char.'); 
grid on;

%Impulzová charakteristika
subplot(2, 2, 2);
impulse(system_P, 'b', system_Z, 'r--', system_P_c2d, 'g:');
title('Impulzová char.'); 
grid on;

%Amplitudová charakteristika
subplot(2, 2, 3);
h = bodeplot(system_P, system_Z, system_P_c2d);
setoptions(h, 'FreqUnits', 'rad/s', 'PhaseVisible', 'off'); 
title('Amplitudová char.'); 
grid on;

%Fázová charakteristika
subplot(2, 2, 4);
h = bodeplot(system_P, system_Z, system_P_c2d);
setoptions(h, 'FreqUnits', 'rad/s', 'MagVisible', 'off'); 
title('Fázová char.'); 
grid on;

%% --- ÚKOL 4: Realizace číslicového filtru (Diferenční rovnice) ---

B_coeffs = [0, num_b1, num_b0];
A_coeffs = [1, den_a1, den_a0];


y_loop = zeros(size(signal)); % Inicializace výstupního pole nulami
x = signal; % Pro přehlednost v rovnici

for k = 3:length(x)
    y_loop(k) = B_coeffs(2)*x(k-1) + B_coeffs(3)*x(k-2) ...
              - A_coeffs(2)*y_loop(k-1) - A_coeffs(3)*y_loop(k-2);
end

y_filter = filter(B_coeffs, A_coeffs, signal);

figure('Name', 'Úkol 1: Realizace cislicoveho filtru');


%původní signál s rušením
subplot(2, 1, 1);
plot(t, signal, 'Color', [0.7 0.7 0.7]);
hold on;
title('Původní zašuměný signál');
xlabel('Čas [s]'); ylabel('Amplituda');
grid on;

subplot(2, 1, 2);

plot(t, y_loop, 'b', 'LineWidth', 1.5); 
hold on;

plot(t, y_filter, 'g--', 'LineWidth', 2);

plot(t, f_filtrovana, 'r', 'LineWidth', 1);

legend('For cyklus', 'Funkce filter()', 'Úkol 1');
title('Porovnání filtrovaných signálů');
xlabel('Čas [s]'); ylabel('Amplituda');
grid on;

%zoom
%xlim([0 0.01]);