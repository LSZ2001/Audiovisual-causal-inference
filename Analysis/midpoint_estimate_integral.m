%% Evaluates a definite integral using midpoint estimate.
% For usage, seeMidpointMethodTest.m.
function [definite_integral] =  midpoint_estimate_integral(fun, a, b, num_bins)
    % Uses midpoint estimate to estimate definite integral of fun from a to b.
    % fun must be an anonymous function of 1 variable only.

    binedges = linspace(a,b,num_bins+1);
    binwidth =  binedges(2)-binedges(1);
    midpoints = binedges(1:end-1) + binwidth/2;
	definite_integral = 0;
    for i=1:length(midpoints)
        definite_integral = definite_integral + binwidth*(fun(midpoints(i)));
    end
end