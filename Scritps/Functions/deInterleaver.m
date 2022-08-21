%==========================================================================
%                      Block Interleaver
%==========================================================================
%   x   --> Secuencia de datos.
%   n   --> Profundidad del interleaver.
%
%   y   --> Secuencia entrelazada (en formato fila).
%==========================================================================
function [y,Block] = deInterleaver(x,n)
    dim = size(x);
    if(dim(1)>1)    % Se revisa que ingrese vector fila.
        x = x.';
    end
    if(mod(length(x),n)~=0)
        x = [x (1:mod(length(x),n))*0];
    end
    Block = reshape(x,[n length(x)/n]);
    y = reshape(Block.',1,[]);  %Se traspone aca para dejar bien formateada la matriz Block.
end