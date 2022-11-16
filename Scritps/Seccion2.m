%==========================================================================
%                            TRABAJO FINAL
%                      Curso Wireless - CoMyS 2022
%                          Llorente, J. F. 
%               Desempeño con CSI parcial y tasa fija...
%==========================================================================
addpath('./Functions');
clc;    clear variables;
% close all;
%============================CONFIGURACION=================================
theta = 0;      REP_CODE_FLAG = 1;  INT_CODE_FLAG = 1;
LW = 2;       ts = 5e-6;           M = 16 ;    
%==========================================================================
N=log2(M);
% NumS=M*250;
% EbN0=15;%[dB]
% A=1; %Es=A^2 y Eb=Es ---> Analíticos
% bits_t=randi([0 1],1,NumS*N); 
% %bk=sign(rand(1,M)-0.5); %Lleva mas tiempo de ejecución.
% [ak,]=generarSimbolos(bits_t,A,M);
% Es=var(ak);
% Eb_1=Es/N;
% N0_veces=Eb_1/(10^(EbN0/10));
% WGNi = sqrt(N0_veces/2)*randn(1,NumS); %RBG con varianza N0/2.
% WGNq = sqrt(N0_veces/2)*randn(1,NumS); % En modelo son ind entonces genero dos veces.
% ak_ni = ak+WGNi;
% ak_nq = WGNq;
% ak_n = cos(theta)*ak_ni-sin(theta)*ak_nq + 1i*(ak_ni*sin(theta)+ak_nq*cos(theta));
% 
% figure; plot(real(ak),0*ak_n,'bo','LineWidth',LW),grid on,set(gca,'FontSize',12);
% title('Constelación BPSK');
% xlabel('$$\phi_0(t)$$','Interpreter','Latex'),xlim(A*[-2 2]);
% xticks(A*[-2 -1 0 1 2]);
% xticklabels({'-2A','-A','0','A','2A'});
% 
% figure; plot(ak_n,'bo','LineWidth',LW),grid on,set(gca,'FontSize',12);
% title(sprintf("Scatter Plot de BPSK con Eb/N0 = %gdB",EbN0));
% xlabel('$$\phi_0(t)$$','Interpreter','Latex'),xlim(A*[-2 2]),ylim(A*[-2 2]);
% xticks(A*[-2 -1 0 1 2]),yticks(A*[-2 -1 0 1 2]);
% xticklabels({'-2A','-A','0','A','2A'}),yticklabels({'-2A','-A','0','A','2A'});

%% Estimación de la PEB
Tc = 0.0162;
samps_inTc = floor(Tc/ts);

NumB=5e7;
Rs = 200e3;
T = 1;
n = 1;
if (REP_CODE_FLAG == 1)
    n = 4;
%     T = n*T;
end
samps_toDecorr = 0*20*samps_inTc + 4;
samps_toDecorr = samps_toDecorr - mod(samps_toDecorr,n);

Ns_xloop = Rs;
Nb_xloop=Ns_xloop*N;
loop = ceil(NumB/Nb_xloop);
NumS = loop*Ns_xloop;
NumB = loop*Nb_xloop;

paso=2;     limite=40; %Parametros para la relevación de la curva.
EsN0_dB=0:paso:limite;
Peb=0.*EsN0_dB;

Es = 1;
A = SymbEnergy2Amp(M,Es);
[Asignacion_bits, Asignacion_coords]=AsignacionBITSyCOORD(M,A); %Por si quiero modificar el código entre bits y simbs.
jj=1;
for EsN0db=0:paso:limite
    p = (1:loop)*0;
    %Señal. Cada loop tiene una realización de bits distinta.
    bits_t = randi([0 1],1,Nb_xloop); %Bits transmitidos.

%     Es = *var(c_noise);
%     A = SymbEnergy2Amp(M,Es);
    
    [aki,akq] = generarSimbolos(bits_t,A,M);
    ak_t = aki + 1i*akq;
    if (REP_CODE_FLAG == 1)
        ak_rep = repCod(ak_t,n);
        ak = ak_rep;    %Para debuggin.
    else
        ak = ak_t;
    end
    
    if (INT_CODE_FLAG == 1)
        [ak_int,ceros] = Interleaver(ak,samps_toDecorr);
        ak = ak_int;    %Para debuggin.
    else
        ceros = 0;
    end
%     A = sqrt(var(ak)/2);
    %Ruido. Se decide generar uno distinto en cada cambio de SNR y no en
    %cada loop porque eleva bastante el costo computacional (tiempo de
    %simulación).
    N0_veces = var(ak)/(10^(EsN0db/10));
    if (INT_CODE_FLAG ~=1)
        WGNi = sqrt(N0_veces/2)*randn(1,Ns_xloop*n); %RBG con varianza N0/2 = 1/2.
        WGNq = sqrt(N0_veces/2)*randn(1,Ns_xloop*n); % En modelo son ind entonces genero dos veces.
        c_noise = (WGNi + 1i*WGNq);
    else
        WGNi = sqrt(N0_veces/2)*randn(1,n*Ns_xloop + ceros); %RBG con varianza N0/2 = 1/2.
        WGNq = sqrt(N0_veces/2)*randn(1,n*Ns_xloop + ceros); % En modelo son ind entonces genero dos veces.
        c_noise = (WGNi + 1i*WGNq);
        
    end
    
    for iteracion=1:loop    
        %Canal. Cada loop tiene una realización de canal distinta.
        if (INT_CODE_FLAG ~= 1)
            h = CanalFlat(n*T,ts);
            h_mat = reshape(h,n,[]);
        else
%             h_mat = zeros(n,floor(T/ts)+1);
            h = CanalFlat(n*T + ceros*ts,ts);
%             h_mat(i,:) = CanalFlat(T,ts);
            h_mat = reshape(deInterleaver(h,samps_toDecorr),n,[]);
            norm_hmat = abs(h_mat).^2;
            norm_hmat = sqrt(sum(norm_hmat));   %A esta altura ya tengo la norma de cada c = [c0 c1 c2 c3] para cada simbolo repetido.
        end
%         h = 0*h + 1;
        
        y_n = ak.*h + c_noise;
        
%         norm_h_aux = reshape(h,n,[]);
%         norm = sqrt(sum(abs(norm_h_aux).^2));
%         norm_h = repCod(norm,n);
        
%         if(EsN0db == limite - paso)
%             debug = 1;
%         end
%     Simbolos_r_1=ReceptorOptimo(real(y_n)./(channel_50segs),imag(y_n)./(channel_50segs),A_1,M,Asignacion_coords_1);
%     Dividiendo por Channel se cancelan los modulos y las fases se restan
%     por lo que deja de estar presente la secuencia de ganancias de canal.
%         h_mat = reshape(h,n,[]);  %Si no estoy generando 4 canales indps.
        if (REP_CODE_FLAG == 1)
            if (INT_CODE_FLAG ==1)
                y_n = deInterleaver(y_n,samps_toDecorr);
            end
            y_mrl = reshape(y_n,n,[]).*conj(h_mat);
            y_mrl = sum(y_mrl)./norm_hmat;
            Simbolos_r_1=ReceptorOptimo(real(y_mrl),imag(y_mrl),A*norm_hmat,M,Asignacion_coords);
        else
            Simbolos_r_1=ReceptorOptimo(real(y_n./h),imag(y_n./h),A,M,Asignacion_coords);
        end
        
        % Esta decisión de pasar amplitudes pesadas por h es para que el
        % ruido solo sufra un cambio de fase y no cambio de amplitud.
%         if (INT_CODE_FLAG == 1)
%             Simbolos_r_1 = deInterleaver(Simbolos_r_1,n);
%         end
%         if (REP_CODE_FLAG == 1) %Es como si los simbolos ya hubiesen sido deInterleaved.
%             Simbolos_r_1 = repDeco(Simbolos_r_1,n,M,Asignacion_coords);    %Simbolos ya decodificados.
%         end
        
        bits_r_1=ConvaBits(Simbolos_r_1,Asignacion_coords,M);
        p(iteracion)=sum(bits_r_1(1:(end-ceros))~=bits_t)/Nb_xloop; 
%         p(iteracion)=sum(Simbolos_r_1 ~= ak_t)/Ns_xloop;
    end
    Peb(jj) = mean(p);
    jj=jj+1;
end
% Pes = Peb*N;
%% Gráficos
fprintf("Simulación Wireless por canal con desvanecimiento Rayleigh.\nEsquema de modulación: %s",ModScheme(M));
if(REP_CODE_FLAG == 1)
    fprintf(" + código de repetición de %d veces",n);
end
if(INT_CODE_FLAG == 1)
    fprintf(" + entrelazado de profundidad %d.\n",n);
else
    fprintf(".\n");
end
fprintf("%g bits simulados.\n",NumB);

EsN0_veces = 10.^(EsN0_dB/10);
Peb_BPSK=qfunc(sqrt(2*EsN0_veces));%Igual a Pes
Peb_QPSK=qfunc(sqrt(2*EsN0_veces/2));  %Dos BPSK ind en fase y cuadratura.
Pes_16QAM_holgada=3*qfunc(sqrt(4/5*EsN0_veces/4));
Peb_16QAM_holgada=Pes_16QAM_holgada/4;
Peb_BPSK_fading = 0.5*(1-sqrt(EsN0_veces./(1+EsN0_veces)));
Peb_QPSK_fading = 0.5*(1-sqrt(EsN0_veces./(2+EsN0_veces)));
Peb_16QAM_fading = 5/2./EsN0_veces;
Pendiente = 1./(4*EsN0_veces).^n;
Peb_BPSK_fading_REPyINT = nchoosek(2*n-1,n).*Pendiente;
figure;
switch M
    case 2
        if(REP_CODE_FLAG && INT_CODE_FLAG)
            semilogy(EsN0_dB,Peb,EsN0_dB,Peb_BPSK,EsN0_dB,Peb_BPSK_fading_REPyINT,'--k','LineWidth',LW/4);
        else
            semilogy(EsN0_dB,Peb,EsN0_dB,Peb_BPSK,EsN0_dB,Peb_BPSK_fading,'--k','LineWidth',LW/4);
        end
        legend('Relevada','Teórica AWGN','Teórica FADING');
    case 4
        if(REP_CODE_FLAG && INT_CODE_FLAG)
            semilogy(EsN0_dB,Peb,EsN0_dB,Peb_QPSK,EsN0_dB,Pendiente,'--k','LineWidth',LW/4);
            legend('Relevada','Teórica AWGN','Pendiente teórica FADING');
        else
            semilogy(EsN0_dB,Peb,EsN0_dB,Peb_QPSK,EsN0_dB,Peb_QPSK_fading,'--k','LineWidth',LW/4);
            legend('Relevada','Teórica AWGN','Teórica FADING');
        end
        
    otherwise
        if(REP_CODE_FLAG && INT_CODE_FLAG)
            semilogy(EsN0_dB,Peb,EsN0_dB,Peb_16QAM_holgada,EsN0_dB,Pendiente,'--k','LineWidth',LW/4);
            legend('Relevada','Teórica AWGN','Pendiente teórica FADING');
        else
            semilogy(EsN0_dB,Peb,EsN0_dB,Peb_16QAM_holgada,EsN0_dB,Peb_16QAM_fading,'--k','LineWidth',LW/4);
            legend('Relevada','Teórica AWGN','Asintótica FADING');
        end
end
set(gca,'FontSize',11);
title(sprintf("Curva de probabilidad de error de bit %s (%g bits)",ModScheme(M),NumB));
grid on, ylabel('BER','Interpreter','Latex'),xlabel('$$E_s/N_0 [dB]$$','Interpreter','Latex');
ylim([9.9e-6 1]);