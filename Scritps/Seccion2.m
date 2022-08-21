%==========================================================================
%                            TRABAJO FINAL
%                      Curso Wireless - CoMyS 2022
%                          Llorente, J. F. 
%               Desempeño con CSI parcial y tasa fija...
%==========================================================================
addpath('./Functions');
clc;    clear variables; close all;
%============================CONFIGURACION=================================
theta = 0;
LW = 2;       ts = 5e-6;           M = 4 ;    
%==========================================================================
N=log2(M);
NumS=M*250;
EbN0=15;%[dB]
A=1; %Es=A^2 y Eb=Es ---> Analíticos
bits_t=randi([0 1],1,NumS*N); 
%bk=sign(rand(1,M)-0.5); %Lleva mas tiempo de ejecución.
[ak,]=generarSimbolos(bits_t,A,M);
Es=var(ak);
Eb_1=Es/N;
N0_veces=Eb_1/(10^(EbN0/10));
WGNi = sqrt(N0_veces/2)*randn(1,NumS); %RBG con varianza N0/2.
WGNq = sqrt(N0_veces/2)*randn(1,NumS); % En modelo son ind entonces genero dos veces.
ak_ni = ak+WGNi;
ak_nq = WGNq;
ak_n = cos(theta)*ak_ni-sin(theta)*ak_nq + 1i*(ak_ni*sin(theta)+ak_nq*cos(theta));

figure; plot(real(ak),0*ak_n,'bo','LineWidth',LW),grid on,set(gca,'FontSize',12);
title('Constelación BPSK');
xlabel('$$\phi_0(t)$$','Interpreter','Latex'),xlim(A*[-2 2]);
xticks(A*[-2 -1 0 1 2]);
xticklabels({'-2A','-A','0','A','2A'});

figure; plot(ak_n,'bo','LineWidth',LW),grid on,set(gca,'FontSize',12);
title(sprintf("Scatter Plot de BPSK con Eb/N0 = %gdB",EbN0));
xlabel('$$\phi_0(t)$$','Interpreter','Latex'),xlim(A*[-2 2]),ylim(A*[-2 2]);
xticks(A*[-2 -1 0 1 2]),yticks(A*[-2 -1 0 1 2]);
xticklabels({'-2A','-A','0','A','2A'}),yticklabels({'-2A','-A','0','A','2A'});

%% Estimación de la PEB
NumB=1e6;
Rs = 200e3;

Ns_xloop = Rs;
Nb_xloop=Ns_xloop*N;
loop = floor(NumB/Nb_xloop);
NumS = loop*Ns_xloop;

paso=2;     limite=40; %Parametros para la relevación de la curva.
EsN0_dB=0:paso:limite;
Peb=0.*EsN0_dB;

jj=1;
for EsN0db=0:paso:limite
    p = (1:loop)*0;
    %Ruido. Se decide generar uno distinto en cada cambio de SNR y no en
    %cada loop porque eleva bastante el costo computacional (tiempo de
    %simulación).
    N0_veces = 1;
    WGNi = sqrt(N0_veces/2)*randn(1,Ns_xloop); %RBG con varianza N0/2 = 1/2.
    WGNq = sqrt(N0_veces/2)*randn(1,Ns_xloop); % En modelo son ind entonces genero dos veces.
    c_noise = (WGNi + 1i*WGNq);
    for iteracion=1:loop    
        %Canal. Cada loop tiene una realización de canal distinta.
        h = CanalFlat(1,ts);
%         h = 0*h + 1;
        %Señal. Cada loop tiene una realización de bits distinta.
        Es = (10^(EsN0db/10))*var(c_noise);
        A = SymbEnergy2Amp(M,Es);
        bits_t=randi([0 1],1,Nb_xloop); %Bits transmitidos.
        [Asignacion_bits, Asignacion_coords]=AsignacionBITSyCOORD(M,A); %Por si quiero modificar el código entre bits y simbs.
        [aki,akq] = generarSimbolos(bits_t,A,M);
        ak = aki + 1i*akq;
        y_n = ak.*h + c_noise;
%     Simbolos_r_1=ReceptorOptimo(real(y_n)./(channel_50segs),imag(y_n)./(channel_50segs),A_1,M,Asignacion_coords_1);
%     Dividiendo por Channel se cancelan los modulos y las fases se restan
%     por lo que deja de estar presente la secuencia de ganancias de canal.
        Simbolos_r_1=ReceptorOptimo(real(y_n).*conj(h)./abs(h),imag(y_n).*conj(h)./abs(h),A,M,Asignacion_coords);
        bits_r_1=ConvaBits(Simbolos_r_1,Asignacion_coords,M);
        p(iteracion)=sum(bits_r_1~=bits_t)/Nb_xloop; 
    end
    Peb(jj) = mean(p);
    jj=jj+1;
end
%% Graficos
fprintf("Sistema %s",ModScheme(M));
    
EsN0_veces = 10.^(EsN0_dB/10);
Peb_BPSK=qfunc(sqrt(2*EsN0_veces));%Igual a Pes
Peb_BPSK_fading = 0.5*(1-sqrt(EsN0_veces./(1+EsN0_veces)));
figure;semilogy(EsN0_dB,Peb,EsN0_dB,Peb_BPSK,EsN0_dB,Peb_BPSK_fading,'--k','LineWidth',LW/4),legend('Relevada','Teórica AWGN','Teórica FADING'),set(gca,'FontSize',11);
title(sprintf("Curva de probabilidad de error de bit %s (%g bits)",ModScheme(M),NumB));
grid on, ylabel('BER','Interpreter','Latex'),xlabel('$$E_s/N_0 [dB]$$','Interpreter','Latex');
ylim([9.9e-6 1]);