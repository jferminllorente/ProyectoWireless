function System = ModScheme(M)
    switch (M)
        case 2
            System = "BPSK";
        case 4
            System = "QPSK";
        case 16
            System = "16QAM";
        otherwise
            System = "Error";
    end
end