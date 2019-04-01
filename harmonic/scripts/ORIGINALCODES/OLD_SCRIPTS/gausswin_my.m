function w = gausswin(L, a)

error(nargchk(1,2,nargin,'struct'));

% Default value for Alpha
if nargin < 2 || isempty(a),
    a = 2.5;
end

% Check for valid window length (i.e., N < 0)
% [L,w,trivialwin] = check_order(L);
% if trivialwin, return, end;

% Compute window according to [1]
N = L-1;
n = (0:N)'-N/2;
w = exp(-(1/2)*(a*n/(N/2)).^2);





