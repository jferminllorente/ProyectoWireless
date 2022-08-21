clc;    clear variables; close all;
LW = 2;
%% Estimación de la PEB
theta16 = 0;
theta = 0;
A_1 = 1;    A_2 = sqrt(1/2);    A_3 = sqrt(1/10);
M_1 = 2;    N_1 = log2(M_1);
M_2 = 4;    N_2 = log2(M_2);
M_3 = 16;   N_3 = log2(M_3);
NumB=1e6;
bits_t=randi([0 1],1,NumB); %Bits transmitidos.

NumS_2=NumB/N_2;
NumS_3=NumB/N_3;

[Asignacion_bits_2, Asignacion_coords_2]=AsignacionBITSyCOORD(M_2,A_2);
[Asignacion_bits_3, Asignacion_coords_3]=AsignacionBITSyCOORD(M_3,A_3);

[ak_2_i,ak_2_q]=generarSimbolos(bits_t,A_2,M_2);
ak_2=ak_2_i+1i*ak_2_q;
Simbolos_t_2=ReceptorOptimo(real(ak_2),imag(ak_2),A_2,M_2,Asignacion_coords_2);
Es_2=var(ak_2);
Eb_2=Es_2/N_2;

[ak_3_i,ak_3_q]=generarSimbolos(bits_t,A_3,M_3);
ak_3=ak_3_i+1i*ak_3_q;
Simbolos_t_3=ReceptorOptimo(real(ak_3),imag(ak_3),A_3,M_3,Asignacion_coords_3);
Es_3=var(ak_3);
Eb_3=Es_3/N_3;

paso=2;     limite=40; %Parametros para la relevación de la curva.
EbN0_dB=0:paso:limite;
Peb2=0.*EbN0_dB;
Peb3=0.*EbN0_dB;

jj=1;
for EbN0db=0:paso:limite
    h_2 = CanalFlat(2.5,5e-6);
    h_3 = CanalFlat(1.25,5e-6);
    h_2 = h_2*0+1;
    h_3 = h_3*0+1;
    N0_veces_2=Eb_2/(10^(EbN0db/10)); %Como elegí el valor de A de manera
%     tal que la Eb=1 siempre puedo usar un solo valor de EbN0.
    N0_veces_3=Eb_3/(10^(EbN0db/10));

    WGNi_2 = sqrt(N0_veces_2/2)*randn(1,NumS_2); 
    WGNi_3 = sqrt(N0_veces_3/2)*randn(1,NumS_3); 
    WGNq_2 = sqrt(N0_veces_2/2)*randn(1,NumS_2); 
    WGNq_3 = sqrt(N0_veces_3/2)*randn(1,NumS_3); 
    c_noise_2 = complex(WGNi_2,WGNq_2);
    c_noise_3 = complex(WGNi_3,WGNq_3);
    
%     ak_2_ni=ak_2_i+WGNi_2;
%     ak_2_nq=ak_2_q+WGNq_2;
    ak_2_n = (ak_2_i + 1i*ak_2_q).*h_2 +  c_noise_2; 
%     ak_3_ni=ak_3_i+WGNi_3;
%     ak_3_nq=ak_3_q+WGNq_3;
    ak_3_n = (ak_3_i + 1i*ak_3_q).*h_3 +  c_noise_3;

    Simbolos_r_2=ReceptorOptimo(real(ak_2_n).*conj(h_2)./abs(h_2),imag(ak_2_n).*conj(h_2)./abs(h_2),A_2,M_2,Asignacion_coords_2);
    Simbolos_r_3=ReceptorOptimo(real(ak_3_n).*conj(h_3)./abs(h_3),imag(ak_3_n).*conj(h_3)./abs(h_3),A_3,M_3,Asignacion_coords_3);
    
    bits_r_2=ConvaBits(Simbolos_r_2,Asignacion_coords_2,M_2);
    bits_r_3=ConvaBits(Simbolos_r_3,Asignacion_coords_3,M_3);
    Peb2(jj)=sum(bits_r_2~=bits_t)/NumB;
    Peb3(jj)=sum(bits_r_3~=bits_t)/NumB;

    jj=jj+1;
end


%% Comparación con valores teóricos
EbN0_veces=10.^(EbN0_dB/10);
Peb_BPSK=qfunc(sqrt(2*EbN0_veces));%Igual a Pes
Peb_QPSK_fading = 0.5*(1-sqrt(EbN0_veces./(2+EbN0_veces)));
% Pes_QPSK=2*Peb_BPSK;
Peb_QPSK=Peb_BPSK;  %Dos BPSK ind en fase y cuadratura.
% Pes_16QAM_holgada=15*qfunc(sqrt(4/5*EbN0_veces));
% Peb_16QAM_holgada=Pes_16QAM_holgada/N_3;
% Pes_QPSK_holgada=3*Peb_BPSK;    %Integro teniendo en cuenta a todos los simbolos con 3 simbolos a 2A. Cota de la unión.
% Peb_QPSK_holgada=Pes_QPSK_holgada/N_2;

figure;semilogy(EbN0_dB,Peb2,EbN0_dB,Peb_BPSK,EbN0_dB,Peb_QPSK_fading,'--k','LineWidth',LW/4),legend('Relevada','Teórica AWGN','Teórica FADING'),set(gca,'FontSize',11);
title(sprintf("Curva de probabilidad de error de bit %s (%g bits)",ModScheme(4),NumB));
grid on, ylabel('BER','Interpreter','Latex'),xlabel('$$E_s/N_0 [dB]$$','Interpreter','Latex');
ylim([9.9e-6 1]);