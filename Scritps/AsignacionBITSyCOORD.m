%==========================================================================
%                     ASIGNACION de SIMBOLOS y COORDENADAS
%                           x Fermín Llorente
%==========================================================================
%
%   M   ---> Cantidad de simbolos del sistema.
%   A   ---> Valor de la amplitud A en la constelación.
%
%  [Asignacion_bits , Asignacion_coords]=AsignacionBITSyCOORD(M,A);
function [Asignacion_bits , Asignacion_coords]=AsignacionBITSyCOORD(M,A)
switch M
    case 2
        s0=0;   coordS0=-A;
        s1=1;   coordS1=A;
        Asignacion_bits=[s0;s1];
        Asignacion_coords=[coordS0,coordS1];
    case 4
        s0=[0 0];   coordS0=A+1i*A;
        s1=[0 1];   coordS1=-A+1i*A;
        s3=[1 1];   coordS3=-A-1i*A;
        s2=[1 0];   coordS2=A-1i*A;
        Asignacion_bits=[s0;s1;s2;s3];
        Asignacion_coords=[coordS0,coordS1,coordS2,coordS3];
    case 16
        s0=[0 0 0 0];   coordS0=-3+31i;
        s1=[0 0 0 1];   coordS1=-3+1i;
        s2=[0 0 1 0];   coordS2=-3-31i;
        s3=[0 0 1 1];   coordS3=-3 -1i;
        s4=[0 1 0 0];   coordS4=-1+31i;
        s5=[0 1 0 1];   coordS5=-1 +1i;
        s6=[0 1 1 0];   coordS6=-1-31i;
        s7=[0 1 1 1];   coordS7=-1-1i;
        s8=[1 0 0 0];   coordS8=3+31i;
        s9=[1 0 0 1];   coordS9=3+1i;
        s10=[1 0 1 0];  coordS10=3-31i;
        s11=[1 0 1 1];  coordS11=3-1i;
        s12=[1 1 0 0];  coordS12=1+31i;
        s13=[1 1 0 1];  coordS13=1+1i;
        s14=[1 1 1 0];  coordS14=1-31i;
        s15=[1 1 1 1];  coordS15=1-1i;
        
        Asignacion_bits=[s0;s1;s2;s3;s4;s5;s6;s7;s8;s9;s10;s11;s12;s13;s14;s15];
        Asignacion_coords=A*[coordS0,coordS1,coordS2,coordS3,coordS4,coordS5,coordS6,coordS7,coordS8,...
                                coordS9,coordS10,coordS11,coordS12,coordS13,coordS14,coordS15];
end