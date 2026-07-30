function boundaries_current = Update_boundaries(boundaries_all, t)
%UPDATE_BOUNDARIES Return currently active boundaries at time t.
%
% boundaries_current = UPDATE_BOUNDARIES(boundaries_all, t) filters the
% struct array boundaries_all to only those boundaries that exist at the 
% given time t.
%
% A boundary is considered to exist at time t if:
% - appear_timestamp <= t AND
% - disappear_timestamp == -1 (always exists) OR disappear_timestamp >= t
%
% Inputs:
%   boundaries_all: Struct array of boundaries, as returned by Load_boundaries
%   t: Current time (numeric)
%
% Output:
%   boundaries_current: Struct array of boundaries active at time t

    % Logical index of boundaries that are active at time t
    activeIdx = [boundaries_all.appear_timestamp] <= t & ( [boundaries_all.disappear_timestamp] >= t | [boundaries_all.disappear_timestamp] == -1 );

    % Filter boundaries_all based on activeIdx
    boundaries_current = boundaries_all(activeIdx);
end
