%==========================================================================
%                          Block Interleaver 2
%                                 JFL
%   Se considera diversidad n ya que es como tener n 
%   canales en tiempo, lo que logra que la pendiente de la curva de 
%   Peb vs. SNR caiga con SNR^n
%==========================================================================
%           y = Interleaver2(x,n)
%   x   --> Secuencia de datos a entrelazar.
%   n   --> Profundidad del interleaver.
%
%   y   --> Secuencia entrelazada (en formato fila).
%==========================================================================
function [y,ceros] = Interleaver2(x,n,corrSamples)
    dim = size(x);
    if(dim(1)>1)    % Se revisa que ingrese vector fila.
        x = x.';
    end
    if(mod(length(x),corrSamples)~=0)
        x = [x ((1:mod(length(x),n))*0)];    %Completo con ceros.
    end
%     x_mat = reshape(x,corrSamples,[]).';    %corrSamples elementos en cada fila tengo.
%     tam = size(x_mat);
%     for i=1:tam(1)
    Block = reshape(x,[corrSamples length(x)/corrSamples]).';
    y = reshape(Block,1,[]);
end