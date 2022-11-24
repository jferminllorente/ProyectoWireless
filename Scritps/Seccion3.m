%==========================================================================
%                            TRABAJO FINAL
%                      Curso Wireless - CoMyS 2022
%                          Llorente, J. F. 
%                    Simulación a tasa variable
%==========================================================================
addpath('./Functions');
clc;    clear variables; close all;
%% %============================CONFIGURACION==============================
LW = 2;       ts = 5e-6;  
INTERVAL_SET = 1;   INTERVAL_CENTER = 2;    INTERVAL_END = 3;
Interval_OP = [1 1 0; 1 2 1; 1 0 -1;0 1 0];    OP = INTERVAL_SET;

NONE  = 0;  %   - No se transmite nada.                    (0)
BPSK4 = 1;  %   - BPSK4 : BPSK con código de repetición 4. (1)  
QPSK4 = 2;  %   - QPSK4 : QPSK con código de repetición 4. (2)    
QPSK2 = 3;  %   - QPSK2 : QPSK con código de repetición 2. (3)
QPSK  = 4;  %   - QPSK  : QPSK sin codigo de repetición.   (4)
QAM16 = 5;  %   - QAM16 : 16QAM sin código de repetición.  (5)
%==========================================================================
%% Estimación de la PEB
NumB=1e7;
Rs = 200e3;
T = 50;

T_c = 0.018;
samples_in_Tc = round(T_c/ts);
loop = floor(NumB/samples_in_Tc);

EsN0dB_vect = 0:40;
p = EsN0dB_vect*0;
R = EsN0dB_vect*0;
for jj = 1:length(EsN0dB_vect)
    times0_NONE = 0;    times1_BPSK4 = 0;    times2_QPSK4 = 0;   %Para debug. 
    times3_QPSK2 = 0;    times4_QPSK = 0;     times5_QAM16 = 0;
    
    h = CanalFlat(2*T,ts);
    EsN0dB = EsN0dB_vect(jj);     
    EsN0veces = 10^(EsN0dB/10);
    bits_t=randi([0 1],1,NumB);
    Bindx = 1;  ii=1;
    bits_r = [];
    while (Bindx<=(NumB-samples_in_Tc*4) && ii<=loop)
        i = ii - times0_NONE;   %Para evitar que no se transmitan los datos a los que les toca un canal pinchado.

        indx_c = floor( (Interval_OP(2,OP)*(ii-Interval_OP(1,OP)) + ...
            Interval_OP(4,OP)) *samples_in_Tc/Interval_OP(2,OP)) + ...
            Interval_OP(3,OP);   %Indice para tomar el valor en el inicio, medio o final del intervalo de largo T_c.
        SNReff = 20*log10(abs(h(indx_c))) + EsN0dB;
        SNRrange = (SNReff<-10)*NONE + (SNReff>=-10 && SNReff<-5)*BPSK4 + ...
            (SNReff>=-5 && SNReff<0)*QPSK4 + (SNReff>=0 && SNReff<5)*QPSK2 + ...
            (SNReff>=5 && SNReff<10)*QPSK + (SNReff>=10)*QAM16 ;
        switch SNRrange
            case NONE  % No se transmite nada si es menor que -10dB.
                times0_NONE = times0_NONE + 1;
                aux_bits = [];
            case BPSK4  % BPSK con código de repetición 4 entre -10 y -5dB.
                [aux_bits,Bindx] = EtEwirelessComm(bits_t,h((ii-1)*samples_in_Tc+1:ii*samples_in_Tc),Bindx,BPSK4,EsN0dB);
                times1_BPSK4 = times1_BPSK4 + 1;
            case QPSK4  % QPSK con código de repetición 4 entre -5 y 0dB.
                times2_QPSK4 = times2_QPSK4 + 1;
                [aux_bits,Bindx] = EtEwirelessComm(bits_t,h((ii-1)*samples_in_Tc+1:ii*samples_in_Tc),Bindx,QPSK4,EsN0dB);
            case QPSK2  % QPSK con código de repetición 2 entre 0 y 5dB.
                times3_QPSK2 = times3_QPSK2 + 1;
                [aux_bits,Bindx] = EtEwirelessComm(bits_t,h((ii-1)*samples_in_Tc+1:ii*samples_in_Tc),Bindx,QPSK2,EsN0dB);
            case QPSK  % QPSK entre 5 y 10dB.
                times4_QPSK = times4_QPSK + 1;
                [aux_bits,Bindx] = EtEwirelessComm(bits_t,h((ii-1)*samples_in_Tc+1:ii*samples_in_Tc),Bindx,QPSK,EsN0dB);
            otherwise % 16-QAM si es mayor que 10dB.
                times5_QAM16 = times5_QAM16 + 1;
                [aux_bits,Bindx] = EtEwirelessComm(bits_t,h((ii-1)*samples_in_Tc+1:ii*samples_in_Tc),Bindx,QAM16,EsN0dB);
        end
        ii = ii + 1;
        bits_r = [bits_r aux_bits];
    end
    % Calculo de tasa media relevada.
    time_sim = (ii-1)*samples_in_Tc*ts; %Ultimo indice tomado del canal, multiplicado por Ts para tener el tiempo.
    R(jj) = length(bits_r)/time_sim;    %Tasa media para este valor de SNR.
    bits_tt = bits_t(1:length(bits_r));  %Se descartan los bits que no se transmitieron...
    p(jj) = sum(bits_r~=bits_tt)/length(bits_tt);

end
EsN0_veces = 10.^(EsN0dB_vect/10);
% Calculo de tasa media teórica.
P_ab = @(a,b,SNR) exp(-a./SNR).*(a./SNR + 1) - exp(-b./SNR).*(b./SNR + 1);
R_m = 1/ts*(1/4*P_ab(10^(-1),10^(-0.5),EsN0_veces) + 1/2*P_ab(10^(-0.5),10^(0),EsN0_veces) ... 
    + P_ab(10^(0),10^(0.5),EsN0_veces) + 2*P_ab(10^(0.5),10^(1),EsN0_veces) + 4*P_ab(10^(1),10^(100),EsN0_veces));

%% Gráficos
figure;
semilogy(EsN0dB_vect,p,'LineWidth',LW/4);
set(gca,'FontSize',11);
% title("Curva de probabilidad de error de bit tasa variable");
grid on, ylabel('BER','Interpreter','Latex'),xlabel('$$E_s/N_0 [dB]$$','Interpreter','Latex');
ylim([9.9e-6 1]);
figure;
plot(EsN0dB_vect,R,EsN0dB_vect,R_m,'--k','LineWidth',LW/4);
set(gca,'FontSize',11);
% title("Tasa media obtenida");
grid on, ylabel('Rate [bits/s]','Interpreter','Latex'),xlabel('$$E_s/N_0 [dB]$$','Interpreter','Latex');
legend('Relevada','Teórica');