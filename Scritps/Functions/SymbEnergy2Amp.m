function A = SymbEnergy2Amp(M,Es)
    switch (M)
        case 2
            A = sqrt(Es);
        case 4
            A = sqrt(Es/2);
        otherwise    %16QAM
            A = sqrt(Es/10);
    end
end