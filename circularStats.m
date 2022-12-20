function [m, s] = circularStats(x,per)
% [m, s] = circularStats(x,per)
%
% Calculate circular statistics (mean and variance) 
%
% Inputs:
%   x: samples
%   per: period
%
% Outputs:
%   m: mean
%   s: standard deviation

n = length(x);
 xRad = x/per*2*pi;
 C = sum(cos(xRad));
 Cbar = C/n;
 S = sum(sin(xRad));
 Sbar = S/n;
 
 R = sqrt(C^2+S^2);
 Rbar = R/n;
 
 T = atan2(Sbar,Cbar);
 if T<0
     T = T+2*pi;
 end
 m = T/2/pi*per;
 
 v = sqrt(-2*log(Rbar));
 s = v/2/pi*per;