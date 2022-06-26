function [result, I, C, DampedNaturalFrequencies] = GetPeakT_partD(c1, c2, I1, I2, FirstDampedDisk, SecondDampedDisk, I, C, K, k, ObservedDisk)
    I(end-1,end-1)=I1;
    I(end,end)=I2;
    Iinv=I^-1;
    
    C(end-1, end-1)=c1;
    C(end, end)=c2;
    C(end-1, FirstDampedDisk)=-c1;
    C(end, SecondDampedDisk)=-c2;
    C(FirstDampedDisk, FirstDampedDisk)=c1;
    C(FirstDampedDisk, end-1)=-c1;
    C(SecondDampedDisk, SecondDampedDisk)=c2;
    C(SecondDampedDisk, end)=-c2;
    
    % Calculate damped natural frequencies:
    A=[zeros(size(I)), eye(size(I)); -Iinv*K, -Iinv*C];
    DampedNaturalFrequencies=eig(A)/1i;
    DampedNaturalFrequencies=real(DampedNaturalFrequencies);
    DampedNaturalFrequencies=DampedNaturalFrequencies(DampedNaturalFrequencies>0);
    
    %Get T for the damped natural frequencies:
    result=zeros(size(DampedNaturalFrequencies));
    for i=1:length(result)
        result(i)=GetT(DampedNaturalFrequencies(i), I, C, K, k, ObservedDisk);
    end
end