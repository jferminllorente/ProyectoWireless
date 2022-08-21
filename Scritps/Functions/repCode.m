%Siempre devuelve vector fila.
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