% This code is the minmod function used in the slope limiter
% zhiwen zhang
% 20130513

function res = minmod(a,b)

%% parameter 
%% compute
res = 0.5*(sign(a) + sign(b)).*min(abs(a),abs(b));