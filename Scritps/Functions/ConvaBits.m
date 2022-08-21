%==========================================================================
%                      Conversion de simbolos a bits
%                           x Fermín Llorente
%==========================================================================
%
%   s   ---> Vector de simbolos 
%   Asignacion_coords   ---> Asignación de coordenadas para cada simbolo.
%   M   ---> Cantidad de simbolos del sistema.
%
%   b = ConvaBits(s,Asignacion_coords,M);
function b=ConvaBits(s,Asignacion_coords,M)
    N=log2(M);
    dim=length(s);
    simbs=zeros(M,dim);

    switch M
        case 2
            b=(s==Asignacion_coords(2));%S1 es un uno y S0 es un cero.
        case 4
%             simbs(1,:)=(s==Asignacion_coords(1))*0;
            simbs(1,:)=s.*0; 
            simbs(2,:)=(s==Asignacion_coords(2))*1;
            simbs(3,:)=(s==Asignacion_coords(3))*2;
            simbs(4,:)=(s==Asignacion_coords(4))*3;
            simbs_num=sum(simbs);
            bits_array=fliplr(de2bi(simbs_num));
            b=reshape(bits_array.',[1 dim*N]);
        case 16
%             simbs(1,:)=(s==Asignacion_coords(1))*0;
            simbs(1,:)=s.*0;
            simbs(2,:)=(s==Asignacion_coords(2))*1;
            simbs(3,:)=(s==Asignacion_coords(3))*2;
            simbs(4,:)=(s==Asignacion_coords(4))*3;
            simbs(5,:)=(s==Asignacion_coords(5))*4;
            simbs(6,:)=(s==Asignacion_coords(6))*5;
            simbs(7,:)=(s==Asignacion_coords(7))*6;
            simbs(8,:)=(s==Asignacion_coords(8))*7;
            simbs(9,:)=(s==Asignacion_coords(9))*8;
            simbs(10,:)=(s==Asignacion_coords(10))*9;
            simbs(11,:)=(s==Asignacion_coords(11))*10;
            simbs(12,:)=(s==Asignacion_coords(12))*11;
            simbs(13,:)=(s==Asignacion_coords(13))*12;
            simbs(14,:)=(s==Asignacion_coords(14))*13;
            simbs(15,:)=(s==Asignacion_coords(15))*14;
            simbs(16,:)=(s==Asignacion_coords(16))*15;
            simbs_num=sum(simbs);
            bits_array=fliplr(de2bi(simbs_num));
            b=reshape(bits_array.',[1 dim*N]);
    end
end