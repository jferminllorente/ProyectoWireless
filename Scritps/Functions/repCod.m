ot%==========================================================================
%           Codificador de código de repetición de orden n               |¦
%                                 JFL                                    |¦
%==========================================================================
%           y = repCod(x,n)                                              |¦
%   x   --> Secuencia de bits sin codificar.                             |¦
%   n   --> Cantidad de veces que se repite cada bit.                    |¦
%                                                                        |¦
%   y   --> Secuencia codificada (en formato fila).                      |¦
%==========================================================================
function y = repCod(x,n)
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