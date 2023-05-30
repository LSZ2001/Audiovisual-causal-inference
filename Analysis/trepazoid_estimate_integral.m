%% Evaluates a definite integral using midpoint estimate.
% For usage, seeMidpointMethodTest.m.
function [definite_integral] =  trepazoid_estimate_integral(fun, a, b, num_bins)
    % Uses midpoint estimate to estimate definite integral of fun from a to b.
    % fun must be an anonymous function of 1 variable only.

    binedges = linspace(a,b,num_bins+1);
    binwidth =  binedges(2)-binedges(1);
	definite_integral = 0;
    for i=1:(length(binedges)-1)
        definite_integral = definite_integral + binwidth.*((fun(binedges(i))+fun(binedges(i+1)))./2);
    end
end