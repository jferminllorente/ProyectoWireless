%==========================================================================
%                          Block Interleaver
%                                 JFL
%   Se elige profundidad n de manera tal que se obtenga diversidad 4 ya que
%   es como tener 4 canales en tiempo, lo que logra que la pendiente de la 
%   curva de Peb vs. SNR caiga con SNR^4.
%==========================================================================
%           y = Interleaver(x,n)
%   x   --> Secuencia de datos a entrelazar.
%   n   --> Profundidad del interleaver. mmm
%
%   y   --> Secuencia entrelazada (en formato fila).
%==========================================================================
function [y,ceros] = Interleaver(x,n)
    dim = size(x);
    ceros = 0;
    if(dim(1)>1)    % Se revisa que ingrese vector fila.
        x = x.';
    end
    if(mod(length(x),n)~=0)
        ceros = (n-mod(length(x),n));
        x = [x (1:(n-mod(length(x),n)))*0];    %Completo con ceros.
    end
    Block = reshape(x,[length(x)/n n]).';
    y = reshape(Block,1,[]);
end