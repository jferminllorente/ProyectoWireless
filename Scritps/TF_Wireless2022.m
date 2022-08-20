clc;    clear variables;    close all;
%==========================================================================
%                            TRABAJO FINAL
%                      Curso Wireless - CoMyS 2022
%                          Llorente, J. F. 
%==========================================================================

%% Introducción
T = 1; %Tiempo de simulación en segundos.
CanalFlat;
vel_ms = V;
D_s = 2*V/lambda; %Doppler spread.
T_c = 1/(4*D_s);    %Tiempo de coherencia.
% Ancho de banda de coherencia se calcula como el soporte de A_C(deltaF) =
% TF(en retardo tau){A_c(tau)} --> Filmina 12 de clase 2...

% Gráfica de la ganancia de canal
figure;
plot((0:N-1)*ts, abs(h))
grid on, grid minor;
title('Ganancia del canal')
ylabel('Valor absoluto (veces)')
xlabel('Tiempo (s)')

% Gráfica del espectro Doppler (puede tardar un rato...)
[DOP,ff]=pwelch(h,[],[],2*fr,fs);
S_a_neg = (ff>-fDm & ff<fDm).*1./(pi*fDm*sqrt(1-(ff/fDm).^2));
% S_a = (ff<fDm)./sqrt(1-(ff/fDm).^2);
figure;
plot(ff,DOP,'k'),grid on, grid minor,hold on;
title('Espectro de la ganancia de canal'),ylabel('Densidad Espectral de potencia [W/Hz]'),xlabel('Frecuencia(Hz)');
ylim([-0.1*max(DOP) max(DOP)*1.1]);
plot(ff,S_a_neg,'--r');
legend('Realización','Teórico');

%% Realizaciones y obtener PDP
L = 100; %Realizaciones.
C = zeros(length(h),L);
C(:,1) = h;
for i=2:L
    CanalFlat;
    C(:,i) = h;
end
A_c = mean(abs(C).^2,2);    %A_c(\tau)
figure;
plot((0:N-1)*ts, A_c)
grid on, grid minor;
title('A_c(\tau)')
ylabel('Valor absoluto (veces)')
xlabel('Tiempo (s)')
%% Primer enfoque
p_t = A_c/(sum(A_c)*ts);
figure;
plot((0:N-1)*ts, p_t)
grid on, grid minor;
title('p_T(\tau)')
ylabel('Valor absoluto (veces)')
xlabel('Tiempo (s)')

A_C_deltaf = fftshift(fft(A_c));
df = fs/length(A_c);
f = -fs/2 : df : fs/2 - df;
figure;
plot(f,abs(A_C_deltaf)),grid on, grid minor;
title('A_C(\Delta f)')
ylabel('Valor absoluto (veces)')
xlabel('Tiempo (s)')

%% Segundo enfoque
A_C_deltaf2_realizations = zeros(length(h)*2-1,L);
Cf = fftshift(fft(C,[],1)*df);
for i = 1:L
    A_C_deltaf2_realizations(:,i) = xcorr(Cf(:,i),Cf(:,i));
end
A_C_deltaf2 = mean(A_C_deltaf2_realizations,2);
[~,deltaf] = xcorr(Cf(:,1),Cf(:,1));
figure;
plot(deltaf,abs(A_C_deltaf2)),grid on, grid minor;
title('A_C(\Delta f) segundo método')
ylabel('Valor absoluto (veces)')
xlabel('Tiempo (s)')

%% Verificar Tiempo de Coherencia de manera empirica.
[AC_empirica,deltat] = xcorr(h,h);
AC_empirica = AC_empirica*ts;
figure;
plot(deltat*ts,abs(AC_empirica)),grid on,grid minor;

% El Tc da 0.00812*2 = 0.0162, si tomamos el Tc = 1/Bc = 0.01 y es aprox
% igual. Es arbitrario elegri el dividido 4 o no, en Goldsmith no tiene
% div 4.