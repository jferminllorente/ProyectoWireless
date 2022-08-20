clc;    clear variables;    close all;

%% FORO
N = 1000;
for snr=1:50
    snrnum=10^(snr/10);
    psk_fading(snr)=0.5*(1-sqrt(snrnum/(1+snrnum)));    %Se obtiene la Peb teorica.
end
snr=1:50;
figure;
semilogy(snr,psk_fading,'y');
grid on;
datanrzpolar=randsrc(1,N);       % Se obtiene un vector de 1 y -1 aleatorio.
for snr=1:50
    signal_amp=sqrt(10^(snr/10))*datanrzpolar;  %Se le asigna la energia de simbolo a cada simbolo segun SNR.
    for mcr=1:100
        noise=randn(1,N);
        noise1=randn(1,N);
        cnoise=complex(noise,noise1);
        cnoise=cnoise/(sqrt(var(cnoise)));  %Se normaliza para que tenga potencia unitaria, de manera que la SNR es igual a Es.
        p_f=2*pi*randn(1,N);
        % p_f = 0*(1:N);
        %=======================================
        %       Se crea el canal con fading
        %=======================================
        x1=randn(1,N);
        x2=randn(1,N);
        x1=x1/(sqrt(var(x1)));
        x2=x2/(sqrt(var(x2)));
        for i=1:N
            alpha(i)=sqrt(x1(i)^2+x2(i)^2);
        end
        chi=alpha.^2;
        chi_mean=mean(chi);
        alpha_normalised=alpha./sqrt(chi_mean);
        for i=1:N
            h(i)=alpha_normalised(i)*complex(cos(p_f(i)),-sin(p_f(i)));
            % h(i)=alpha_normalised(i);
            r(i)=h(i)*signal_amp(i)+cnoise(i);
        end
        for i=1:N
            mrc(i)=r(i)*conj(h(i));
        end
        for i=1:N
            if(real(mrc(i))>=0)
                decision(i)=1;
            else
                decision(i)=-1;
            end
        end
        hamm=0;
        for i=1:N
            if(decision(i)~=datanrzpolar(i))
                hamm=hamm+1;
            end
        end
        pe(mcr)=hamm/N;
    end
    pepsk(snr)=mean(pe);
end
snr=1:50;
hold on;
semilogy(snr,pepsk,'r*');
grid on;

%% Test SNR 
