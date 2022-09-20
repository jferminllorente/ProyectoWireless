%==========================================================================
%                         Block Deinterleaver
%                                 JFL
%==========================================================================
%           y = deInterleaver(x,n)
%   x   --> Secuencia de datos entrelazada.
%   n   --> Profundidad con la que se entrelazaron los datos.
%
%   y   --> Secuencia desentrelazada (en formato fila).
%==========================================================================
function y = deInterleaver(x,n)
    dim = size(x);
    if(dim(1)>1)    % Se revisa que ingrese vector fila.
        x = x.';
    end
    Block = reshape(x,[n length(x)/n]);
    y = reshape(Block.',1,[]);  %Se traspone aca para dejar bien formateada la matriz Block.
end