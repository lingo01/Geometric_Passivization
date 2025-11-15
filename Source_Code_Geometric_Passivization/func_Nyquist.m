function [Hvals, ReVals, ImVals, omega] = func_Nyquist(G, Nsample)
%FUNC_NYQUIST_FIXED Compute frequency response and always return omega vector
%   [H, Re, Im, omega] = func_Nyquist_fixed(G, omega)
%   If omega omitted, defaults to linspace(0,100,2001).

omega = 10.^linspace(-3, 3, Nsample);

% Ensure omega is row
omega = reshape(omega, 1, []);

% Compute response
H = squeeze(freqresp(G, omega));
H = reshape(H, 1, []);

Hvals = H;
ReVals = real(Hvals);
ImVals = imag(Hvals);
end
