function [ dist ] = periodic_dist_calc(dimes,sysLen,selectId,pConf,pSelect )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
dimes; % to shut debugger up
delx = (pConf-pSelect);
delx = abs(delx);
delx_compli = abs(delx - sysLen) ;
delx = min(delx_compli,delx);
delx = delx';
dist = sum(delx.^2).^0.5;
if selectId~=0
    dist(selectId) = [] ;
end
end

