function [ rn ] = trunc_gauss_sample( center, sig, lb, ub )
%trunc_gauss_sample Sample from truncated normal distribution
%   Detailed explanation goes here

notFound = 1;
rn = nan(length(center),1);
mask = logical(zeros(length(center),1));

while notFound
    % Draw where still needed
    draw = randn(sum(~mask),1)*sig + center(~mask);
    % Gather in return array
    rn(~mask) = draw;
    % Check whether return array already complies to boundaries
    mask = (rn>=lb) & (rn<=ub);
    % One single violation requires a redraw on the entries where needed
    if mask
        notFound = 0;
    end
end

end

