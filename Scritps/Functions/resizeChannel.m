%==========================================================================
%                         Desechar muestras canal
%                                 JFL
%==========================================================================
%   h   --> Canal de largo N.
%   num --> Cantidad de muestras que quiero de 'h'.
%   OP  --> Punto donde quiero que esten centradas las muestras de mi
%   canal nuevo.
% 
%   reh --> Canal con num muestras de N, centradas en OP.
%==========================================================================
function reh = resizeChannel(h,num,OP)
INTERVAL_SET = 1;   INTERVAL_CENTER = 2;    INTERVAL_END = 3;
    switch OP
        case INTERVAL_SET
            reh = h(1:num);
        case INTERVAL_CENTER
            if mod(num,2)
                reh = h(floor(length(h)/2)-floor(num/2) : floor(length(h)/2) + floor(num/2));
            else
                reh = h(floor(length(h)/2)-floor(num/2) : floor(length(h)/2) + floor(num/2) - 1);
            end
        case INTERVAL_END
            reh = h((length(h)-num):end);
    end
end