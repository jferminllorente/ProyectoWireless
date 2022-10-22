%==========================================================================
%                     End to End Wireless Communication
%                                 JFL
%==========================================================================
%           y = EEwirelessComm(x,c,SCHEME)
%   bits_t   --> Secuencia de bits a enviar.
%   h        --> Secuencia de amplitudes del canal.
%   SCHEME   --> Esquema de modulación + codificación.
%   - BPSK4 : BPSK con código de repetición 4. (1)
%   - QPSK4 : QPSK con código de repetición 4. (2)    
%   - QPSK2 : QPSK con código de repetición 2. (3)
%   - QPSK  : QPSK sin codigo de repetición.   (4)
%   - QAM16 : 16QAM sin código de repetición.  (5)
%   EsN0dB   --> SNR en dB de los simbolos transmitidos.
%   Peb        --> Secuencia entrelazada (en formato fila).
%==========================================================================
function [bits_r,Bindx_out] = EtEwirelessComm(bits_t,h,Bindx,SCHEME,EsN0dB,OP)
Es = 1;
Interval_OP = [1 1 0; 1 2 1; 1 0 -1;0 1 0];
switch SCHEME
    case 1
        repN = 4;
        M = 2;  N = log2(M);
        A = SymbEnergy2Amp(M,Es);
        [~, Asignacion_coords]=AsignacionBITSyCOORD(M,A);
        NumS = length(h);
        bits_actuales = bits_t(Bindx:(Bindx + floor(NumS*N/repN) - 1));
        if (mod(length(bits_actuales),N))
            bits_actuales = bits_actuales(1:end-mod(length(bits_actuales),N));
        end
        Bindx_out = Bindx + length(bits_actuales);
        [aki,akq] = generarSimbolos(bits_actuales,A,M);
        ak = aki + 1i*akq;
        ak = repCod(ak,repN);
        if length(ak) ~= length(h)
            h = h(1:length(ak));    %Matcheo largos si NumS/N/rep no fuese un entero...
        end
        NumS = length(ak);
        N0_veces = var(ak)/(10^(EsN0dB/10));
        WGNi = sqrt(N0_veces/2)*randn(1,NumS); %RBG con varianza N0/2 = 1/2.
        WGNq = sqrt(N0_veces/2)*randn(1,NumS); % En modelo son ind entonces genero dos veces.
        c_noise = (WGNi + 1i*WGNq);
        y_n = ak.*h + c_noise;
        Simbolos_r = ReceptorOptimo(real(y_n./h),imag(y_n./h),A,M,Asignacion_coords);
        bits_r = ConvaBits(Simbolos_r,Asignacion_coords,M);
%         Peb = sum(bits_r_1~=bits_t)/NumB;
    case 2
%         indx_c = floor( (Interval_OP(2,OP)*(ii-Interval_OP(1,OP)) + Interval_OP(4,OP)) *samples_in_Tc/Interval_OP(2,OP)) + Interval_OP(3,OP);
%         h_trunc = resizeChannel(h,length(ak),OP);   %Ya tengo al canal del mismo largo que los simbolos a transmitir.
    case 3
        
    case 4
        repN = 1;
        M = 4;  N = log2(M);
        A = SymbEnergy2Amp(M,Es);
        [~, Asignacion_coords]=AsignacionBITSyCOORD(M,A);
        NumS = length(h);
        bits_actuales = bits_t(Bindx:(Bindx + floor(NumS*N/repN) - 1));
        if (mod(length(bits_actuales),N))
            bits_actuales = bits_actuales(1:end-mod(length(bits_actuales),N));
        end
        Bindx_out = Bindx + length(bits_actuales);
        [aki,akq] = generarSimbolos(bits_actuales,A,M);
        ak = aki + 1i*akq;
%         ak = repCod(ak,repN);
        if length(ak) ~= length(h)
            h = h(1:length(ak));    %Matcheo largos si NumS/N/rep no fuese un entero...
        end
        NumS = length(ak);
        N0_veces = var(ak)/(10^(EsN0dB/10));
        WGNi = sqrt(N0_veces/2)*randn(1,NumS); %RBG con varianza N0/2 = 1/2.
        WGNq = sqrt(N0_veces/2)*randn(1,NumS); % En modelo son ind entonces genero dos veces.
        c_noise = (WGNi + 1i*WGNq);
        y_n = ak.*h + c_noise;
        Simbolos_r = ReceptorOptimo(real(y_n./h),imag(y_n./h),A,M,Asignacion_coords);
        bits_r = ConvaBits(Simbolos_r,Asignacion_coords,M);
    case 5
        repN = 1;
        M = 16;  N = log2(M);
        A = SymbEnergy2Amp(M,Es);
        [~, Asignacion_coords]=AsignacionBITSyCOORD(M,A);
        NumS = length(h);
        bits_actuales = bits_t(Bindx:(Bindx + floor(NumS*N/repN) - 1));
        if (mod(length(bits_actuales),N))
            bits_actuales = bits_actuales(1:end-mod(length(bits_actuales),N));
        end
        Bindx_out = Bindx + length(bits_actuales);
        [aki,akq] = generarSimbolos(bits_actuales,A,M);
        ak = aki + 1i*akq;
%         ak = repCod(ak,repN);
        if length(ak) ~= length(h)
            h = h(1:length(ak));    %Matcheo largos si NumS/N/rep no fuese un entero...
        end
        NumS = length(ak);
        N0_veces = var(ak)/(10^(EsN0dB/10));
        WGNi = sqrt(N0_veces/2)*randn(1,NumS); %RBG con varianza N0/2 = 1/2.
        WGNq = sqrt(N0_veces/2)*randn(1,NumS); % En modelo son ind entonces genero dos veces.
        c_noise = (WGNi + 1i*WGNq);
        y_n = ak.*h + c_noise;
        Simbolos_r = ReceptorOptimo(real(y_n./h),imag(y_n./h),A,M,Asignacion_coords);
        bits_r = ConvaBits(Simbolos_r,Asignacion_coords,M);
    otherwise 
        disp('El esquema esta mal!');
end
end