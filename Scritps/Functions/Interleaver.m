%==========================================================================
%                          Block Interleaver
%                                 JFL
%   Se considera diversidad n (con n la profundidad) ya que es como tener n 
%   canales en tiempo, lo que logra que la pendiente de la curva de 
%   Peb vs. SNR caiga con SNR^n
%==========================================================================
%           y = Interleaver(x,n)
%   x   --> Secuencia de datos a entrelazar.
%   n   --> Profundidad del interleaver.
%
%   y   --> Secuencia entrelazada (en formato fila).
%==========================================================================
function y = Interleaver(x,n)
    dim = size(x);
    if(dim(1)>1)    % Se revisa que ingrese vector fila.
        x = x.';
    end
    if(mod(length(x),n)~=0)
        x = [x ((1:mod(length(x),n))*0+x(end))];    %Completo con el ultimo simbolo.
    end
    Block = reshape(x,[length(x)/n n]).';
    y = reshape(Block,1,[]);
end