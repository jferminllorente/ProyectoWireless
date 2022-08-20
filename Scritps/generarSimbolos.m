%==========================================================================
%                           Generar simbolos
%                           x Fermín Llorente
%==========================================================================
%   Función que recibe los bits a transmitir y le asigna coordenadas a cada
%   combinación de bits. Devuelve las coordenadas complejas de cada simbolo
%   transmitido.
%
%   b   ---> Vector de bits
%   A   ---> Valor de la amplitud A en la constelación.
%   M   ---> Cantidad de simbolos del sistema de comunicaciones.
%
%   [aki,akq]=generarSimbolos(b,A,M)
function [aki,akq]=generarSimbolos(b,A,M)
N=log2(M);
NumS=length(b)/N;
ak=zeros(1,NumS);
switch M
    case 2  %BPSK
        s0=0;   coordS0=-A;
        s1=1;   coordS1=A;
        simbs_0=b==s0;
        simbs_1=b==s1;
        ak=simbs_0*coordS0+simbs_1*coordS1;
    case 4  %QPSK
        % Asigno secuencia de bits y coordenadas a simbolos. 
        s0=[0 0];   coordS0=A+1i*A;
        s1=[0 1];   coordS1=-A+1i*A;
        s2=[1 1];   coordS2=-A-1i*A;
        s3=[1 0];   coordS3=A-1i*A;

        simb_bits=reshape(b,[N NumS]);
        aux0=simb_bits(1:N,:)==s0.';
        aux1=simb_bits(1:N,:)==s1.';
        aux2=simb_bits(1:N,:)==s2.';
        aux3=simb_bits(1:N,:)==s3.';
        simbs_0=aux0(1,:)==1 & aux0(2,:)==1;
        simbs_1=aux1(1,:)==1 & aux1(2,:)==1;
        simbs_2=aux2(1,:)==1 & aux2(2,:)==1;
        simbs_3=aux3(1,:)==1 & aux3(2,:)==1;
        ak=simbs_0*coordS0+simbs_1*coordS1+simbs_2*coordS2+simbs_3*coordS3;
 
    case 16  %16QAM
        s0=[0 0 0 0];   coordS0=-3+3*1i;
        s1=[0 1 0 0];   coordS1=-1+3*1i;
        s2=[1 1 0 0];   coordS2=1+3*1i;
        s3=[1 0 0 0];   coordS3=3+3*1i;
        s4=[0 0 0 1];   coordS4=-3+1i;
        s5=[0 1 0 1];   coordS5=-1 +1i;
        s6=[1 1 0 1];   coordS6=1+1i;
        s7=[1 0 0 1];   coordS7=3+1i;
        s8=[0 0 1 1];   coordS8=-3 -1i;
        s9=[0 1 1 1];   coordS9=-1-1i;
        s10=[1 1 1 1];  coordS10=1-1i;
        s11=[1 0 1 1];  coordS11=3-1i;
        s12=[0 0 1 0];  coordS12=-3-3*1i;
        s13=[0 1 1 0];  coordS13=-1-3*1i;
        s14=[1 1 1 0];  coordS14=1-3*1i;
        s15=[1 0 1 0];  coordS15=3-3*1i;
        simb_bits=reshape(b,[N NumS]);
        aux0=simb_bits(1:N,:)==s0.';
        aux1=simb_bits(1:N,:)==s1.';
        aux2=simb_bits(1:N,:)==s2.';
        aux3=simb_bits(1:N,:)==s3.';
        aux4=simb_bits(1:N,:)==s4.';
        aux5=simb_bits(1:N,:)==s5.';
        aux6=simb_bits(1:N,:)==s6.';
        aux7=simb_bits(1:N,:)==s7.';
        aux8=simb_bits(1:N,:)==s8.';
        aux9=simb_bits(1:N,:)==s9.';
        aux10=simb_bits(1:N,:)==s10.';
        aux11=simb_bits(1:N,:)==s11.';
        aux12=simb_bits(1:N,:)==s12.';
        aux13=simb_bits(1:N,:)==s13.';
        aux14=simb_bits(1:N,:)==s14.';
        aux15=simb_bits(1:N,:)==s15.';
        %Algoritmo que sirve para cualquier codigo, aún si no es el de Gray
        simbs_0=aux0(1,:)==1 & aux0(2,:)==1 & aux0(3,:)==1 & aux0(4,:)==1;
        simbs_1=aux1(1,:)==1 & aux1(2,:)==1 & aux1(3,:)==1 & aux1(4,:)==1;
        simbs_2=aux2(1,:)==1 & aux2(2,:)==1 & aux2(3,:)==1 & aux2(4,:)==1;
        simbs_3=aux3(1,:)==1 & aux3(2,:)==1 & aux3(3,:)==1 & aux3(4,:)==1;
        simbs_4=aux4(1,:)==1 & aux4(2,:)==1 & aux4(3,:)==1 & aux4(4,:)==1;
        simbs_5=aux5(1,:)==1 & aux5(2,:)==1 & aux5(3,:)==1 & aux5(4,:)==1;
        simbs_6=aux6(1,:)==1 & aux6(2,:)==1 & aux6(3,:)==1 & aux6(4,:)==1;
        simbs_7=aux7(1,:)==1 & aux7(2,:)==1 & aux7(3,:)==1 & aux7(4,:)==1;
        simbs_8=aux8(1,:)==1 & aux8(2,:)==1 & aux8(3,:)==1 & aux8(4,:)==1;
        simbs_9=aux9(1,:)==1 & aux9(2,:)==1 & aux9(3,:)==1 & aux9(4,:)==1;
        simbs_10=aux10(1,:)==1 & aux10(2,:)==1 & aux10(3,:)==1 & aux10(4,:)==1;
        simbs_11=aux11(1,:)==1 & aux11(2,:)==1 & aux11(3,:)==1 & aux11(4,:)==1;
        simbs_12=aux12(1,:)==1 & aux12(2,:)==1 & aux12(3,:)==1 & aux12(4,:)==1;
        simbs_13=aux13(1,:)==1 & aux13(2,:)==1 & aux13(3,:)==1 & aux13(4,:)==1;
        simbs_14=aux14(1,:)==1 & aux14(2,:)==1 & aux14(3,:)==1 & aux14(4,:)==1;
        simbs_15=aux15(1,:)==1 & aux15(2,:)==1 & aux15(3,:)==1 & aux15(4,:)==1;
        ak = simbs_0*coordS0 + simbs_1*coordS1 + simbs_2*coordS2 + simbs_3*coordS3 + simbs_4*coordS4 + simbs_5*coordS5 +...
             simbs_6*coordS6 + simbs_7*coordS7 + simbs_8*coordS8 + simbs_9*coordS9 + simbs_10*coordS10 + simbs_11*coordS11 +...
             simbs_12*coordS12 + simbs_13*coordS13 + simbs_14*coordS14 + simbs_15*coordS15;
end
aki=A*real(ak);
akq=A*imag(ak);
end