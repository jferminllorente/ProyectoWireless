%==========================================================================
%                            TRABAJO FINAL
%                      Curso Wireless - CoMyS 2022
%                          Llorente, J. F. 
%                    Simulación a tasa variable...
%==========================================================================
addpath('./Functions');
clc;    clear variables; close all;
%% 
%============================CONFIGURACION=================================
LW = 2;       ts = 5e-6;  
INTERVAL_SET = 1;   INTERVAL_CENTER = 2;    INTERVAL_END = 3;
Interval_OP = [1 1 0; 1 2 1; 1 0 -1;0 1 0];    OP = INTERVAL_SET;

NONE  = 0;  %   - No se transmite nada.                    (0)
BPSK4 = 1;  %   - BPSK4 : BPSK con código de repetición 4. (1)  
QPSK4 = 2;  %   - QPSK4 : QPSK con código de repetición 4. (2)    
QPSK2 = 3;  %   - QPSK2 : QPSK con código de repetición 2. (3)
QPSK  = 4;  %   - QPSK  : QPSK sin codigo de repetición.   (4)
QAM16 = 5;  %   - QAM16 : 16QAM sin código de repetición.  (5)

%% Estimación de la PEB
NumB=1e7;
Rs = 200e3;
T = 50;

T_c = 0.018;
% T_c = 0.05;
samples_in_Tc = round(T_c/ts);
loop = floor(NumB/samples_in_Tc);

% contador = 0;
EsN0dB_vect = 0:40;
p = EsN0dB_vect*0;
R = EsN0dB_vect*0;
for jj = 1:length(EsN0dB_vect)
    times0_NONE = 0;    times1_BPSK4 = 0;    times2_QPSK4 = 0;    
    times3_QPSK2 = 0;    times4_QPSK = 0;     times5_QAM16 = 0;
    
    h = CanalFlat(2*T,ts);
    EsN0dB = EsN0dB_vect(jj);     
    EsN0veces = 10^(EsN0dB/10);
    bits_t=randi([0 1],1,NumB);
    Bindx = 1;  ii=1;
    bits_r = [];
    while (Bindx<=(NumB-samples_in_Tc*4) && ii<=loop)
        i = ii - times0_NONE;

        indx_c = floor( (Interval_OP(2,OP)*(ii-Interval_OP(1,OP)) + Interval_OP(4,OP)) *samples_in_Tc/Interval_OP(2,OP)) + Interval_OP(3,OP);   %Indice para tomar el valor en el inicio, medio o final del intervalo de largo T_c.
        SNReff = 20*log10(abs(h(indx_c))) + EsN0dB;
        SNRrange = (SNReff<-10)*NONE + (SNReff>=-10 && SNReff<-5)*BPSK4 + ...
            (SNReff>=-5 && SNReff<0)*QPSK4 + (SNReff>=0 && SNReff<5)*QPSK2 + ...
            (SNReff>=5 && SNReff<10)*QPSK + (SNReff>=10)*QAM16 ;
%         SNRrange = QAM16;
        switch SNRrange
            case NONE  % No se transmite nada si es menor que -10dB.
                % Aca tengo que dividir en vectores de 3600 y transmitir ahi.
                % Quizas con un reshape...
                times0_NONE = times0_NONE + 1;
                aux_bits = [];
            case BPSK4  % BPSK con código de repetición 4 entre -10 y -5dB.
    %             aux_bits = EtEwirelessComm(bits_t((i-1)*samples_in_Tc+1:i*samples_in_Tc),h((i-1)*samples_in_Tc+1:i*samples_in_Tc),BPSK4,EsN0dB,OP);
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
    time_sim = (ii-1)*samples_in_Tc*ts; %Ultimo indice tomado del canal, multiplicado por Ts para tener el tiempo.
    R(jj) = length(bits_r)/time_sim;    %Tasa media para este valor de SNR.
    bits_tt = bits_t(1:length(bits_r));  %Se descartan los bits que no se transmitieron...
    p(jj) = sum(bits_r~=bits_tt)/length(bits_tt);
%     clc;
%     fprintf("Veces que se transmitió con cada sistema:\n");
%     fprintf("Pinchado:    %d\n",times0_NONE);
%     fprintf("BPSK(rep 4): %d\n",times1_BPSK4);
%     fprintf("QPSK(rep 4): %d\n",times2_QPSK4);
%     fprintf("QPSK(rep 2): %d\n",times3_QPSK2);
%     fprintf("QPSK:        %d\n",times4_QPSK);
%     fprintf("16QAM:       %d\n",times5_QAM16);
end

% fprintf("Veces que se transmitió con cada sistema:\n");
% fprintf("Pinchado:    %d\n",times0_NONE);
% fprintf("BPSK(rep 4): %d\n",times1_BPSK4);
% fprintf("QPSK(rep 4): %d\n",times2_QPSK4);
% fprintf("QPSK(rep 2): %d\n",times3_QPSK2);
% fprintf("QPSK:        %d\n",times4_QPSK);
% fprintf("16QAM:       %d\n",times5_QAM16);

% n = 1;
% if (REP_CODE_FLAG == 1)
%     n = 4;
% %     T = n*T;
% end
% Ns_xloop = Rs;
% Nb_xloop=Ns_xloop*N;
% loop = floor(NumB/Nb_xloop);
% NumS = loop*Ns_xloop;
% 
% paso=2;     limite=40; %Parametros para la relevación de la curva.
% EsN0_dB=0:paso:limite;
% Peb=0.*EsN0_dB;
% 
% Es = 1;
% A = SymbEnergy2Amp(M,Es);
% [Asignacion_bits, Asignacion_coords]=AsignacionBITSyCOORD(M,A); %Por si quiero modificar el código entre bits y simbs.
% jj=1;
% for EsN0db=0:paso:limite
%     p = (1:loop)*0;
%     %Señal. Cada loop tiene una realización de bits distinta.
%     bits_t=randi([0 1],1,Nb_xloop); %Bits transmitidos.
% 
% %     Es = *var(c_noise);
% %     A = SymbEnergy2Amp(M,Es);
%     
%     [aki,akq] = generarSimbolos(bits_t,A,M);
%     ak_t = aki + 1i*akq;
%     if (REP_CODE_FLAG == 1)
%         ak_rep = repCod(ak_t,n);
%         ak = ak_rep;    %Para debuggin.
%     else
%         ak = ak_t;
%     end
%     
% %     if (INT_CODE_FLAG == 1)
% %         ak_int = Interleaver(ak,4);
% %         ak = ak_int;    %Para debuggin.
% %     end
% %     A = sqrt(var(ak)/2);
%     %Ruido. Se decide generar uno distinto en cada cambio de SNR y no en
%     %cada loop porque eleva bastante el costo computacional (tiempo de
%     %simulación).
%     N0_veces = var(ak)/(10^(EsN0db/10));
%     if (INT_CODE_FLAG ~=1)
%         WGNi = sqrt(N0_veces/2)*randn(1,Ns_xloop*n); %RBG con varianza N0/2 = 1/2.
%         WGNq = sqrt(N0_veces/2)*randn(1,Ns_xloop*n); % En modelo son ind entonces genero dos veces.
%         c_noise = (WGNi + 1i*WGNq);
%     else
%         noise_mat = zeros(n,Ns_xloop);
%         for i=1:n
%             WGNi = sqrt(N0_veces/2)*randn(1,Ns_xloop); %RBG con varianza N0/2 = 1/2.
%             WGNq = sqrt(N0_veces/2)*randn(1,Ns_xloop); % En modelo son ind entonces genero dos veces.
%             noise_mat(i,:) = (WGNi + 1i*WGNq);
%         end
%         norm_cnoise = abs(noise_mat).^2;
%         norm_cnoise = sqrt(sum(norm_cnoise));   %A esta altura ya tengo la norma de cada c = [c0 c1 c2 c3] para cada simbolo repetido.
% %         norm_cnoise = repmat(norm_cnoise,4,1);
% %         norm_cnoise = reshape(norm_cnoise,1,[]);    %A esta altura tengo la norma repetida.
%         c_noise = reshape(noise_mat,1,[]);
%     end
%     
%     for iteracion=1:loop    
%         %Canal. Cada loop tiene una realización de canal distinta.
%         if (INT_CODE_FLAG ~= 1)
%             h = CanalFlat(n*T,ts);
%             h_mat = reshape(h,n,[]);
%         else
%             h_mat = zeros(n,floor(T/ts)+1);
%             for i=1:n
%                 h_mat(i,:) = CanalFlat(T,ts);
%             end
%             h = reshape(h_mat,1,[]);
%         end
%         norm_hmat = abs(h_mat).^2;
%         norm_hmat = sqrt(sum(norm_hmat));
%         
% %         h = 0*h + 1;
%         
%         y_n = ak.*h + c_noise;
%         
% %         norm_h_aux = reshape(h,n,[]);
% %         norm = sqrt(sum(abs(norm_h_aux).^2));
% %         norm_h = repCod(norm,n);
%         
% %         if(EsN0db == limite - paso)
% %             debug = 1;
% %         end
% %     Simbolos_r_1=ReceptorOptimo(real(y_n)./(channel_50segs),imag(y_n)./(channel_50segs),A_1,M,Asignacion_coords_1);
% %     Dividiendo por Channel se cancelan los modulos y las fases se restan
% %     por lo que deja de estar presente la secuencia de ganancias de canal.
% %         h_mat = reshape(h,n,[]);  %Si no estoy generando 4 canales indps.
%         if (REP_CODE_FLAG == 1)
%             y_mrl = reshape(y_n,n,[]).*conj(h_mat);
%             y_mrl = sum(y_mrl)./norm_hmat;
%             Simbolos_r_1=ReceptorOptimo(real(y_mrl),imag(y_mrl),A*norm_hmat,M,Asignacion_coords);
%         else
%             Simbolos_r_1=ReceptorOptimo(real(y_n./h),imag(y_n./h),A,M,Asignacion_coords);
%         end
%         
%         % Esta decisión de pasar amplitudes pesadas por h es para que el
%         % ruido solo sufra un cambio de fase y no cambio de amplitud.
% %         if (INT_CODE_FLAG == 1)
% %             Simbolos_r_1 = deInterleaver(Simbolos_r_1,n);
% %         end
% %         if (REP_CODE_FLAG == 1) %Es como si los simbolos ya hubiesen sido deInterleaved.
% %             Simbolos_r_1 = repDeco(Simbolos_r_1,n,M,Asignacion_coords);    %Simbolos ya decodificados.
% %         end
%         
%         bits_r_1=ConvaBits(Simbolos_r_1,Asignacion_coords,M);
%         p(iteracion)=sum(bits_r_1~=bits_t)/Nb_xloop; 
% %         p(iteracion)=sum(Simbolos_r_1 ~= ak_t)/Ns_xloop;
%     end
%     Peb(jj) = mean(p);
%     jj=jj+1;
% end
% % Pes = Peb*N;
% %% Gráficos
% fprintf("Simulación Wireless por canal con desvanecimiento Rayleigh.\nEsquema de modulación: %s",ModScheme(M));
% if(REP_CODE_FLAG == 1)
%     fprintf(" + código de repetición de %d veces",n);
% end
% if(INT_CODE_FLAG == 1)
%     fprintf(" + entrelazado de profundidad %d.\n",n);
% else
%     fprintf(".\n");
% end
% fprintf("%g bits simulados.\n",NumB);
% 
% EsN0_veces = 10.^(EsN0_dB/10);
% Peb_BPSK=qfunc(sqrt(2*EsN0_veces));%Igual a Pes
% Peb_QPSK=qfunc(sqrt(2*EsN0_veces/2));  %Dos BPSK ind en fase y cuadratura.
% Pes_16QAM_holgada=3*qfunc(sqrt(4/5*EsN0_veces/4));
% Peb_16QAM_holgada=Pes_16QAM_holgada/4;
% Peb_BPSK_fading = 0.5*(1-sqrt(EsN0_veces./(1+EsN0_veces)));
% Peb_QPSK_fading = 0.5*(1-sqrt(EsN0_veces./(2+EsN0_veces)));
% Peb_16QAM_fading = 5/2./EsN0_veces;
% Pendiente = 1./(4*EsN0_veces).^n;
% Peb_BPSK_fading_REPyINT = nchoosek(2*n-1,n)./Pendiente;
% figure;
% switch M
%     case 2
%         if(REP_CODE_FLAG && INT_CODE_FLAG)
%             semilogy(EsN0_dB,Peb,EsN0_dB,Peb_BPSK,EsN0_dB,Peb_BPSK_fading_REPyINT,'--k','LineWidth',LW/4);
%         else
%             semilogy(EsN0_dB,Peb,EsN0_dB,Peb_BPSK,EsN0_dB,Peb_BPSK_fading,'--k','LineWidth',LW/4);
%         end
%         legend('Relevada','Teórica AWGN','Teórica FADING');
%     case 4
%         if(REP_CODE_FLAG && INT_CODE_FLAG)
%             semilogy(EsN0_dB,Peb,EsN0_dB,Peb_QPSK,EsN0_dB,Pendiente,'--k','LineWidth',LW/4);
%             legend('Relevada','Teórica AWGN','Pendiente teórica FADING');
%         else
%             semilogy(EsN0_dB,Peb,EsN0_dB,Peb_QPSK,EsN0_dB,Peb_QPSK_fading,'--k','LineWidth',LW/4);
%             legend('Relevada','Teórica AWGN','Teórica FADING');
%         end
%         
%     otherwise
%         if(REP_CODE_FLAG && INT_CODE_FLAG)
%             semilogy(EsN0_dB,Peb,EsN0_dB,Peb_16QAM_holgada,EsN0_dB,Pendiente,'--k','LineWidth',LW/4);
%             legend('Relevada','Teórica AWGN','Pendiente teórica FADING');
%         else
%             semilogy(EsN0_dB,Peb,EsN0_dB,Peb_16QAM_holgada,EsN0_dB,Peb_16QAM_fading,'--k','LineWidth',LW/4);
%             legend('Relevada','Teórica AWGN','Asintótica FADING');
%         end
% end
figure;
semilogy(EsN0dB_vect,p,'LineWidth',LW/4);
set(gca,'FontSize',11);
title("Curva de probabilidad de error de bit tasa variable");
grid on, ylabel('BER','Interpreter','Latex'),xlabel('$$E_s/N_0 [dB]$$','Interpreter','Latex');
ylim([9.9e-6 1]);
figure;
plot(EsN0dB_vect,R,'LineWidth',LW/4);
set(gca,'FontSize',11);
title("Tasa media obtenida");
grid on, ylabel('Rate [bits/s]','Interpreter','Latex'),xlabel('$$E_s/N_0 [dB]$$','Interpreter','Latex');