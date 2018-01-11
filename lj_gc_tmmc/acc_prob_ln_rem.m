function [ lnAccProb ] = acc_prob_ln_rem( T,zz,V,N,dE )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
lnAccProb= -dE/T -log(zz) - log(V) + log(N);
if lnAccProb>0
    lnAccProb = 0;
end

end