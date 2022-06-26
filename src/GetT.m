function result = GetT(omega, I, C, K, k, ObservedDisk)
    alpha=(K-omega^2*I+1i*omega*C)^-1;
    result=abs(alpha(ObservedDisk,1)*k);
end