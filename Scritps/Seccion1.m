clc;    clear variables;    close all;
%==========================================================================
%                            TRABAJO FINAL
%                      Curso Wireless - CoMyS 2022
%                          Llorente, J. F.
%                   Caracterización del modelo de canal.
%==========================================================================

%% Introducción
T = 1; %Tiempo de simulación en segundos.
ts = 5e-6;
h = CanalFlat(T,ts);
fc=900e6;               % Frecuencia de portadora (Hz)
V=60/3.6;  
vel_ms = V;
lambda=3e8/fc;          % Longitud de onda
fDm=V/lambda;           % Doppler máximo
N=round(T/ts);          % # muestras (divisible por 4D)
fs=1/ts;                % Frecuencia de muestreo (Hz)
D=round(fDm*N/fs);      % # de puntos para filtro Jakes (el total es 4D)
fD=D*fs/N;              % Doppler máximo redondeado 
Df=fD/D;               % Paso de frecuencia para la generación del espectro
fr=-fD+Df:Df:fD-Df;    % Eje de frecuencias para filtro Jakes

D_s = 2*V/lambda; %Doppler spread.
T_c = 1/(4*D_s);    %Tiempo de coherencia.

%% Gráfica de la ganancia de canal
figure;
plot((0:N-1)*ts, abs(h))
grid on, grid minor;
title('Ganancia del canal')
ylabel('Valor absoluto (veces)')
xlabel('Tiempo (s)')

%% Gráfica del espectro Doppler (puede tardar un rato...)
[DOP,ff]=pwelch(h,[],[],2*fr,fs);
S_a_neg = (ff>-fDm & ff<fDm).*1./(pi*fDm*sqrt(1-(ff/fDm).^2));
figure;
plot(ff,DOP,'k'),grid on, grid minor,hold on;
title('Espectro de la ganancia de canal'),ylabel('Densidad Espectral de potencia [W/Hz]'),xlabel('Frecuencia(Hz)');
ylim([-0.1*max(DOP) max(DOP)*1.1]);
plot(ff,S_a_neg,'--r');
legend('Realización','Teórico');

%% Verificar Tiempo de Coherencia de manera empirica.
[AC_empirica,deltat] = xcorr(h,h);
AC_empirica = AC_empirica*ts;
figure;
plot(deltat*ts,abs(AC_empirica)),grid on,grid minor;

% El Tc da 0.00812, si tomamos el Tc = 1/Bc = 0.01 y es aprox
% igual. Es arbitrario elegri el dividido 4 o no, en Goldsmith no tiene
% div 4.