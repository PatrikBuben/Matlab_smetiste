clc
clearvars
[signal, frekvence_vzorkovani] = audioread('/Users/patrikbuben/Desktop/BPC-SAS/270566.wav');
audiowrite('/Users/patrikbuben/Desktop/BPC-SAS/270566.wav',signal, frekvence_vzorkovani);
%sound(signal, frekvence_vzorkovani)

%zakladni parametry jako delka a frekvence
N = length(signal); %celkovy pocet vzorku
t = (0:N-1)/frekvence_vzorkovani; %casova osa od 0 do delky signalu [vzorec T = 1/f proste]

%transformace
F = fft(signal); %fourier fast transformation
f = (0:N-1) * (frekvence_vzorkovani / N);

F_abs = abs(F);
f_plot = f(1:floor(N/2));
F_magnituda_plot = F_abs(1:floor(N/2));

figure;
subplot(2, 1, 1)
plot(t, signal)
title('Puvodni signal (ptacci)')
xlabel('Cas [s]')
ylabel('Amplituda')
grid
hold on

subplot(2, 1, 2)
loglog(f_plot(2: end), F_magnituda_plot(2:end));
title('Puvodni amplitudove spektrum (obraz signalu)');
xlabel('Frekvence [Hz]')
ylabel('|F(f)|')
grid
hold on

F_filtrovana = F; %zkopiroval jsem si spektrum, abych ho mohl filtrovat
%ted podle grafu najdeme rusive frekvence ja jsem nasel (peaks kde to
%vystreli jak kokot)
% 10852 a 907.635, 10185.6, 728.567
f_peak1 = 10852; %[Hz]
f_peak2 = 907.635;
f_peak3 = 10185.6;
f_peak4 = 728.567;
f_peak5 = 10900;
f_peak6 = 12301;
f_peak7 = 9162.74;

f_peaks = [f_peak1, f_peak2, f_peak3, f_peak4, f_peak5, f_peak6, f_peak7];
sirka_pasma = 50; %timto odrezavame

for peak_freq = f_peaks

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
loglog(f_plot(2:end), F_filtrovana_abs_plot(2:end), "r");
legend('Puvodni spektrum', 'Filtrovane spektrum');
hold off;

%Tip 3
max_val = max(abs(f_filtrovana));
f_normal = f_filtrovana / max_val;

audiowrite('filtrovano.wav', f_normal, frekvence_vzorkovani);
[novy_signal, nova_fs] = audioread('/Users/patrikbuben/Desktop/BPC-SAS/cvikahome/filtrovano.wav');
sound(novy_signal, nova_fs);