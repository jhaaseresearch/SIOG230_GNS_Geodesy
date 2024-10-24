function [MP1,MP2] = compute_multipath(f1,f2, L1, L2, P1, P2)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
alpha = (f1./f2).^2;
MP1 = P1 - (2/(alpha - 1) + 1).* L1 + (2./(alpha-1)).*L2;
MP2 = P2 - ((2.*alpha)./(alpha - 1)).* L1 + (2./(alpha-1) -1 ).*L2;
end