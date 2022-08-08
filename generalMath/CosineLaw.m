function [gamma] = CosineLaw(a,b,c)
%% calculates the '/gama' angel ("across c edge) using the cosine law
% a,b - edges near the angel
% c - edge across the angel
% retunrs "gama" in radians
% dbstop 8;
%%
m = max([a,b,c]);
a = a./m;
b = b./m;
c = c./m;

cosGamma = (a.*a + b.*b - c.*c) ./ (2.*a.*b);
gamma = acos(cosGamma);
end % function CosineLaw