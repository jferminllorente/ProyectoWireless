%==========================================================================
%                      Block Interleaver
%==========================================================================
%   x   --> Secuencia de datos.
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
        x = [x (1:mod(length(x),n))*0];
    end
    Block = reshape(x,[length(x)/n n]).';
    y = reshape(Block,1,[]);
end