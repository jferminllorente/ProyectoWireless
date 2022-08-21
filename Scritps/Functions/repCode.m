%==========================================================================
%                   Código de repetición de orden n                      |¦
%                                 JFL                                    |¦
%==========================================================================
%           y = repCode(x,n)                                             |¦
%   x   --> Secuencia de datos sin codificar.                            |¦
%   n   --> Cantidad de veces que se repite cada bit.                    |¦
%                                                                        |¦
%   y   --> Secuencia codificada (en formato fila).                      |¦
%==========================================================================
function y = repCode(x,n)
    dim = size(x);
    if(dim(1)>1)
        x = x.';
    end
    aux = zeros(n,length(x));
    for i=1:n
        aux(i,:) = x;
    end
    y = reshape(aux,[1 n*length(x)]);
end