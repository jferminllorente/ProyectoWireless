%==========================================================================
%                            TRABAJO FINAL
%                      Curso Wireless - CoMyS 2022
%                          Llorente, J. F. 
%               Desempeño con CSI parcial y tasa fija...
%==========================================================================
clc;    clear variables; close all;
LW = 2;
theta = 0;
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

figure; plot(real(ak_n),0*ak_n,'bo','LineWidth',LW),grid on,set(gca,'FontSize',12);
title(sprintf("Scatter Plot de BPSK con Eb/N0 = %gdB",EbN0));
xlabel('$$\phi_0(t)$$','Interpreter','Latex'),xlim(A_1*[-2 2]);
xticks(A_1*[-2 -1 0 1 2]);
xticklabels({'-2A','-A','0','A','2A'});

%% Estimación de la PEB
NumB=1e7;
N_1 = 1;
Rs = 200e3;
bits_t=randi([0 1],1,NumB); %Bits transmitidos.
NumS_1=NumB/N_1;
[Asignacion_bits_1, Asignacion_coords_1]=AsignacionBITSyCOORD(M_1,A_1); %Por si quiero modificar el código entre bits y simbs.
[aki_1,~] = generarSimbolos(bits_t,A_1,M_1);
Simbolos_t_1=ReceptorOptimo(real(aki_1),imag(aki_1),A_1,M_1,Asignacion_coords_1); %Simbolos transmitidos.
Es_1=var(aki_1);
Eb_1=Es_1/N_1; %Obtengo la energía de bit real en base a la realización de bits que tengo.
paso=5;     limite=75; %Parametros para la relevación de la curva.
EbN0_dB=0:paso:limite;
Peb1=0.*EbN0_dB;
jj=1;
% T = 50; %Tiempo de simulación para que entren 1e7 bits. 1e7(bits)/200e3(bits/s) = 50 s
channel_50segs = CanalFlat(50,ts);
for EbN0db=0:paso:limite
%     EbN0db = EsN0db - 10*log10(N_1);
    N0_veces_1=Eb_1/(10^(EbN0db/10));

    WGNi_1 = sqrt(N0_veces_1/2)*randn(1,NumS_1); %RBG con varianza N0/2.
    WGNq_1 = sqrt(N0_veces_1/2)*randn(1,NumS_1); % En modelo son ind entonces genero dos veces.
    c_noise = WGNi_1 + 1i*WGNq_1;
    ak_1 = aki_1 + 1i*0;
%     ak_1_ni = ak_1 + WGNi_1;
%     ak_1_nq = WGNq_1;
%     ak_1_n=(cos(theta)*ak_1_ni-sin(theta)*ak_1_nq + 1i*(ak_1_ni*sin(theta)+ak_1_nq*cos(theta)));
    y_n = ak_1.*channel_50segs + c_noise;
%     Simbolos_r_1=ReceptorOptimo(real(y_n)./(channel_50segs),imag(y_n)./(channel_50segs),A_1,M_1,Asignacion_coords_1);
    Simbolos_r_1=ReceptorOptimo(real(y_n).*conj(channel_50segs),imag(y_n).*conj(channel_50segs),A_1,M_1,Asignacion_coords_1);
%     Dividiendo por Channel se cancelan los modulos y las fases se restan
%     por lo que deja de estar presente la secuencia de ganancias de canal.
    bits_r_1=ConvaBits(Simbolos_r_1,Asignacion_coords_1,M_1);
    Peb1(jj)=sum(bits_r_1~=bits_t)/NumB;
    jj=jj+1;
end

%% Graficos
EbN0_veces = 10.^(EbN0_dB/10);
% EbN0_veces = EsN0_veces/N_1;
Peb_BPSK=qfunc(sqrt(2*EbN0_veces));%Igual a Pes
Peb_BPSK_fading = 0.5*(1-sqrt(EbN0_veces./(1+EbN0_veces)));
figure;semilogy(EbN0_dB,Peb1,EbN0_dB,Peb_BPSK,EbN0_dB,Peb_BPSK_fading,'--k','LineWidth',LW/4),legend('Relevada','Teórica AWGN','Teórica FADING'),set(gca,'FontSize',11);
title(sprintf("Curva de probabilidad de error de bit BPSK (%g bits)",NumB));
grid on, ylabel('BER','Interpreter','Latex'),xlabel('$$E_s/N_0 [dB]$$','Interpreter','Latex');
ylim([9.9e-6 1]);


