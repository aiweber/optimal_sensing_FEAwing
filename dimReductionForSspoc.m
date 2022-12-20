 function [w_t, Psi] = dimReductionForSspoc(X, G, SspocPrs)
% [w_t, Psi] = dimReductionForSspoc(X, G, SspocPrs)
% 
% This function reduces the dimensionality of the original data and uses
% LDA to find a decision vector.
%
% Inputs:
%   X: data matrix; nSensorLocs x nDataPts
%   G: vector indicating classes of columns in X; 1 x nDataPts
%   SspocPrs: structure containing parameters for optimization.  Fields:
%     rmodes: number of modes to keep during dim reduction
%     wTrunc: desired number of sensors
%
% Outputs:
%   w_t: decision vector to reconstruct; wTrunc x numClasses-1
%   Psi: reduced data matrix; nSensorLocs x SspocPrs.rmodes
%
% This code relies on the LDA_n function.

[U, Sigma, ~] = svd(X, 0);      % U: nSensorLocs x nSensorLocs; Sigma: nSensorLocs x (n_out*length(rots)); V: (n_out*length(rots)) x (n_out*length(rots))

Psi = U(:, 1:SspocPrs.wTrunc);  % truncate eigenvectors to desired number for dim reduction
a = Psi'*X;                     % tranform X into PC coordinate space (PCs of location); a: rmodes x (n_out*length(rots))
w_t = LDA_n(a, G);                % find decision boundary w (code works with >2D data) 
    
for i = 1:size(w_t,2)                               % normalize truncated w_t vector
    w_t(:,i) = w_t(:,i) / sqrt(w_t(:,i)'*w_t(:,i));
end
