function [ lnAccProb ] = acc_prob_ln_disp(T,dE)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
lnAccProb = -dE/T;

if lnAccProb>0
    lnAccProb = 0;
end

end

