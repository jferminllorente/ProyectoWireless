%==========================================================================
%                            Receptor Óptimo
%                           x Fermín Llorente
%==========================================================================
%
%   Devuelve un vector de largo 'cantidad de simbolos' con las coordenadas
%   de los simbolos que se recibieron luego de tomar la decisión con las
%   fronteras óptimas.
%   
%   ak_r_i ---> Coordenadas recibidas en la rama en fase.
%   ak_r_q ---> Coordenadas recibidas en la rama en cuadratura.
%   A      ---> Amplitud A de las coordenadas de la constelación.
%   M      ---> Numero de simbolos del sistema que se modela.
%   AsignacionCoords ---> Coordenadas de los M simbolos.
%
%   simb=ReceptorOptimo(ak_r_i,ak_r_q,A,M,AsignacionCoords)
function simb=ReceptorOptimo(ak_r_i,ak_r_q,A,M,AsignacionCoords)
    if(imag(ak_r_q)~=0)
        ak_r_q = imag(ak_r_q);
    end
    simb_aux=zeros(M,length(ak_r_i));
    switch M            
        case 2
            simb_aux(1,:)=(ak_r_i<0).*AsignacionCoords(1);
            simb_aux(2,:)=(ak_r_i>0).*AsignacionCoords(2);
            simb=sum(simb_aux);
        case 4
            simb_aux(1,:)=((ak_r_i>0)&(ak_r_q>0)).*AsignacionCoords(1);
            simb_aux(2,:)=((ak_r_i<0)&(ak_r_q>0)).*AsignacionCoords(2);
            simb_aux(3,:)=((ak_r_i>0)&(ak_r_q<0)).*AsignacionCoords(3);
            simb_aux(4,:)=((ak_r_i<0)&(ak_r_q<0)).*AsignacionCoords(4);
            simb=sum(simb_aux);
        case 16 
            simb_aux(1,:)=((ak_r_i<-2*A)&(ak_r_q>2*A)).*AsignacionCoords(1);
            simb_aux(2,:)=((ak_r_i<-2*A)&(ak_r_q>0 & ak_r_q<2*A)).*AsignacionCoords(2);
            simb_aux(3,:)=((ak_r_i<-2*A)&(ak_r_q<-2*A)).*AsignacionCoords(3);
            simb_aux(4,:)=((ak_r_i<-2*A)&(ak_r_q<0 & ak_r_q>-2*A)).*AsignacionCoords(4);
            simb_aux(5,:)=((ak_r_i<0 & ak_r_i>-2*A)&(ak_r_q>2*A)).*AsignacionCoords(5);
            simb_aux(6,:)=((ak_r_i<0 & ak_r_i>-2*A)&(ak_r_q>0 & ak_r_q<2*A)).*AsignacionCoords(6);
            simb_aux(7,:)=((ak_r_i<0 & ak_r_i>-2*A)&(ak_r_q<-2*A)).*AsignacionCoords(7);
            simb_aux(8,:)=((ak_r_i<0 & ak_r_i>-2*A)&(ak_r_q<0 & ak_r_q>-2*A)).*AsignacionCoords(8);
            simb_aux(9,:)=((ak_r_i>2*A)&(ak_r_q>2*A)).*AsignacionCoords(9);
            simb_aux(10,:)=((ak_r_i>2*A)&(ak_r_q>0 & ak_r_q<2*A)).*AsignacionCoords(10);
            simb_aux(11,:)=((ak_r_i>2*A)&(ak_r_q<-2*A)).*AsignacionCoords(11);
            simb_aux(12,:)=((ak_r_i>2*A)&(ak_r_q<0 & ak_r_q>-2*A)).*AsignacionCoords(12);
            simb_aux(13,:)=((ak_r_i>0 & ak_r_i<2*A)&(ak_r_q>2*A)).*AsignacionCoords(13);
            simb_aux(14,:)=((ak_r_i>0 & ak_r_i<2*A)&(ak_r_q>0 & ak_r_q<2*A)).*AsignacionCoords(14);
            simb_aux(15,:)=((ak_r_i>0 & ak_r_i<2*A)&(ak_r_q<-2*A)).*AsignacionCoords(15);
            simb_aux(16,:)=((ak_r_i>0 & ak_r_i<2*A)&(ak_r_q<0 & ak_r_q>-2*A)).*AsignacionCoords(16);
            simb=sum(simb_aux);
    end
end