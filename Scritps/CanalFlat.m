% Simulador de canal con Rayleigh Flat Fading
% Septiembre de 2012 - PAR
% ========================================================================
% h = CanalFlat(T,ts)
% T   --- Tiempo de simulación.
% ts  --- Tiempo de muestreo.
function h = CanalFlat(T,ts)
    fc=900e6;               % Frecuencia de portadora (Hz)
    V=60/3.6;               % velocidad en m/s
    % ========================================================================
    lambda=3e8/fc;          % Longitud de onda
    fDm=V/lambda;           % Doppler máximo
    N=round(T/ts);          % # muestras (divisible por 4D)
    fs=1/ts;                % Frecuencia de muestreo (Hz)
    D=round(fDm*N/fs);      % # de puntos para filtro Jakes (el total es 4D)
    fD=D*fs/N;              % Doppler máximo redondeado 
    % =======================================================================
    Df=fD/D;               % Paso de frecuencia para la generación del espectro
    fr=-fD+Df:Df:fD-Df;    % Eje de frecuencias para filtro Jakes
    jk=sqrt(1./sqrt(1-(fr/fDm).^2));% Filtro en frecuencia!
    JK=fftshift([jk 0]);   % Filtro acomodado para usar FFT

    HC=sqrt(D)*(randn(1,2*D)+1i*randn(1,2*D)); % Espectro gaussiano BLANCO
    HJ2=HC.*JK;             % Espectro con forma JAKES
    HJN=[HC(1:D).*JK(1:D) zeros(1,N-2*D) HC(D+1:2*D).*JK(D+1:2*D)];
    %Completamos la parte nula del espectro 
    %(luego al hacer IFFT se interpolan los valores de ganancia del canal)
    hj=sqrt(N)*ifft(HJN);   % Ganancia del canal compleja CON CORRELACIÓN JAKES
    h=hj./std(hj);          % Normalizamos para tener ganancia media 1
end

%% =======================================================================
% Gráfica de la ganancia de canal
% =======================================================================
% figure;
% plot((0:N-1)*ts, abs(h))
% grid on, grid minor;
% title('Ganancia del canal')
% ylabel('Valor absoluto (veces)')
% xlabel('Tiempo (s)')

%% ======================================================================
% Gráfica del espectro Doppler (puede tardar un rato...)
% =======================================================================

% [DOP,ff]=pwelch(h,[],[],2*fr,fs);
% S_a_neg = (ff>-fDm & ff<fDm).*1./(pi*fDm*sqrt(1-(ff/fDm).^2));
% % S_a = (ff<fDm)./sqrt(1-(ff/fDm).^2);
% figure;
% plot(ff,DOP,'k'),grid on, grid minor,hold on;
% title('Espectro de la ganancia de canal'),ylabel('Densidad Espectral de potencia [W/Hz]'),xlabel('Frecuencia(Hz)');
% ylim([-0.1*max(DOP) max(DOP)*1.1]);
% plot(ff,S_a_neg,'--r');
% legend('Realización','Teórico');

%% Prueba escala logarítmica
% figure;
% plot(ff,10*log10(DOP),'k'),grid on, grid minor,hold on;
% title('Espectro de la ganancia de canal (en escala logarítmica)'),ylabel('Densidad Espectral de potencia [dBw/Hz]'),xlabel('Frecuencia(Hz)');
% % ylim([-0.1*max(DOP) max(DOP)*1.1]);
% plot(ff,10*log10(S_a_neg),'--r');
% legend('Realización','Teórico');