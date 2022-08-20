%==========================================================================
%                            TRABAJO FINAL
%                      Curso Wireless - CoMyS 2022
%                          Llorente, J. F. 
%               Desempeño con CSI parcial y tasa fija...
%==========================================================================
clc;    clear variables; close all;
theta = 0;
LW = 2;
ts = 5e-6;
M_1=2;    N_1=log2(M_1);
NumS=M_1*250;
EbN0=15;%[dB]
A_1=1; %Es=A^2 y Eb=Es ---> Analíticos
bits_t=randi([0 1],1,NumS*N_1); 
%bk=sign(rand(1,M)-0.5); %Lleva mas tiempo de ejecución.
[ak,]=generarSimbolos(bits_t,A_1,M_1);
Es_1=var(ak);
Eb_1=Es_1/N_1;
N0_veces=Eb_1/(10^(EbN0/10));
WGNi = sqrt(N0_veces/2)*randn(1,NumS); %RBG con varianza N0/2.
WGNq = sqrt(N0_veces/2)*randn(1,NumS); % En modelo son ind entonces genero dos veces.
ak_ni = ak+WGNi;
ak_nq = WGNq;
ak_n = cos(theta)*ak_ni-sin(theta)*ak_nq + 1i*(ak_ni*sin(theta)+ak_nq*cos(theta));

figure; plot(real(ak),0*ak_n,'bo','LineWidth',LW),grid on,set(gca,'FontSize',12);
title('Constelación BPSK');
xlabel('$$\phi_0(t)$$','Interpreter','Latex'),xlim(A_1*[-2 2]);
xticks(A_1*[-2 -1 0 1 2]);
xticklabels({'-2A','-A','0','A','2A'});

figure; plot(ak_n,'bo','LineWidth',LW),grid on,set(gca,'FontSize',12);
title(sprintf("Scatter Plot de BPSK con Eb/N0 = %gdB",EbN0));
xlabel('$$\phi_0(t)$$','Interpreter','Latex'),xlim(A_1*[-2 2]),ylim(A_1*[-2 2]);
xticks(A_1*[-2 -1 0 1 2]),yticks(A_1*[-2 -1 0 1 2]);
xticklabels({'-2A','-A','0','A','2A'}),yticklabels({'-2A','-A','0','A','2A'});

%% Estimación de la PEB
NumB=1e7;
N_1 = 1;
Rs = 200e3;

Nb_xloop = Rs;
loop = floor(NumB/Nb_xloop);
NumB = loop*Rs;
Ns_xloop=Nb_xloop/N_1;

paso=2;     limite=40; %Parametros para la relevación de la curva.
EsN0_dB=0:paso:limite;
Peb1=0.*EsN0_dB;

jj=1;
for EsN0db=0:paso:limite
    p = (1:loop)*0;
    %Ruido. Se decide generar uno distinto en cada cambio de SNR y no en
    %cada loop porque eleva bastante el costo computacional (tiempo de
    %simulación).
    N0_veces_1 = 1;
    WGNi_1 = sqrt(N0_veces_1/2)*randn(1,Ns_xloop); %RBG con varianza N0/2 = 1/2.
    WGNq_1 = sqrt(N0_veces_1/2)*randn(1,Ns_xloop); % En modelo son ind entonces genero dos veces.
    c_noise = (WGNi_1 + 1i*WGNq_1);
    for iteracion=1:loop    
        %Canal. Cada loop tiene una realización de canal distinta.
        h = CanalFlat(1,ts);
        %Señal. Cada loop tiene una realización de bits distinta.
        Es_1 = (10^(EsN0db/10))*var(c_noise);
        bits_t=randi([0 1],1,Nb_xloop); %Bits transmitidos.
        [Asignacion_bits_1, Asignacion_coords_1]=AsignacionBITSyCOORD(M_1,sqrt(Es_1)); %Por si quiero modificar el código entre bits y simbs.
        [aki_1,~] = generarSimbolos(bits_t,sqrt(Es_1),M_1);
        ak_1 = aki_1 + 1i*0;
        y_n = ak_1.*h + c_noise;
%     Simbolos_r_1=ReceptorOptimo(real(y_n)./(channel_50segs),imag(y_n)./(channel_50segs),A_1,M_1,Asignacion_coords_1);
%     Dividiendo por Channel se cancelan los modulos y las fases se restan
%     por lo que deja de estar presente la secuencia de ganancias de canal.
        Simbolos_r_1=ReceptorOptimo(real(y_n).*conj(h),imag(y_n).*conj(h),sqrt(Es_1),M_1,Asignacion_coords_1);
        bits_r_1=ConvaBits(Simbolos_r_1,Asignacion_coords_1,M_1);
        p(iteracion)=sum(bits_r_1~=bits_t)/Nb_xloop; 
    end
    Peb1(jj) = mean(p);
    jj=jj+1;
end
%% Graficos
EsN0_veces = 10.^(EsN0_dB/10);
Peb_BPSK=qfunc(sqrt(2*EsN0_veces));%Igual a Pes
Peb_BPSK_fading = 0.5*(1-sqrt(EsN0_veces./(1+EsN0_veces)));
figure;semilogy(EsN0_dB,Peb1,EsN0_dB,Peb_BPSK,EsN0_dB,Peb_BPSK_fading,'--k','LineWidth',LW/4),legend('Relevada','Teórica AWGN','Teórica FADING'),set(gca,'FontSize',11);
title(sprintf("Curva de probabilidad de error de bit BPSK (%g bits)",NumB));
grid on, ylabel('BER','Interpreter','Latex'),xlabel('$$E_s/N_0 [dB]$$','Interpreter','Latex');
ylim([9.9e-6 1]);