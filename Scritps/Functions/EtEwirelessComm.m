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
function [bits_r,Bindx_out] = EtEwirelessComm(bits_t,h,Bindx,SCHEME,EsN0dB)
Es = 1;

BPSK4 = 1;  %   - BPSK4 : BPSK con código de repetición 4. (1)  
QPSK4 = 2;  %   - QPSK4 : QPSK con código de repetición 4. (2)    
QPSK2 = 3;  %   - QPSK2 : QPSK con código de repetición 2. (3)
QPSK  = 4;  %   - QPSK  : QPSK sin codigo de repetición.   (4)
QAM16 = 5;  %   - QAM16 : 16QAM sin código de repetición.  (5)
switch SCHEME
    case BPSK4
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
        
        h_mat = reshape(h,repN,[]); % Se puede ya que h fue llevado al largo de ak que es multiplo de repN.
        norm_hmat = abs(h_mat).^2;
        norm_hmat = sqrt(sum(norm_hmat));
        
        y_mrl = reshape(y_n,repN,[]).*conj(h_mat);
        y_mrl = sum(y_mrl)./norm_hmat;
        
        Simbolos_r = ReceptorOptimo(real(y_mrl),imag(y_mrl),A*norm_hmat,M,Asignacion_coords);
        bits_r = ConvaBits(Simbolos_r,Asignacion_coords,M);

    case QPSK4
        repN = 4;
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
        
        h_mat = reshape(h,repN,[]); % Se puede ya que h fue llevado al largo de ak que es multiplo de repN.
        norm_hmat = abs(h_mat).^2;
        norm_hmat = sqrt(sum(norm_hmat));
        
        y_mrl = reshape(y_n,repN,[]).*conj(h_mat);
        y_mrl = sum(y_mrl)./norm_hmat;
        
        Simbolos_r = ReceptorOptimo(real(y_mrl),imag(y_mrl),A*norm_hmat,M,Asignacion_coords);
        bits_r = ConvaBits(Simbolos_r,Asignacion_coords,M);
    case QPSK2
        repN = 2;
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
        
        h_mat = reshape(h,repN,[]); % Se puede ya que h fue llevado al largo de ak que es multiplo de repN.
        norm_hmat = abs(h_mat).^2;
        norm_hmat = sqrt(sum(norm_hmat));
        
        y_mrl = reshape(y_n,repN,[]).*conj(h_mat);
        y_mrl = sum(y_mrl)./norm_hmat;
        
        Simbolos_r = ReceptorOptimo(real(y_mrl),imag(y_mrl),A*norm_hmat,M,Asignacion_coords);
        bits_r = ConvaBits(Simbolos_r,Asignacion_coords,M);
    case QPSK
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
    case QAM16
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