function [s_hats_PM, tEnd] = midpoint_trepazoid_getpostmean(sigma_fun, x_obs_vals, sigma0,k,sigma_s,mu, a, b, numbins, numerical_method, scale, prior_type, other_spread, w)
    if(nargin==10)
        scale=1; prior_type="SingleGaussian"; otherspread=0; w=0;
    elseif(nargin==11)
        prior_type="SingleGaussian"; otherspread=0; w=0;
    elseif(nargin==12)
        otherspread=0; w=0;
    elseif(nargin==13)
        w=0;
    end
    
    tStart = tic;
    switch prior_type
        case "SingleGaussian"
            prior = @(s) exp(-((s-mu).^2)./(2.*(sigma_s.^2)));
        case "TwoGaussiansOneFixedZero"
        case "TwoGaussiansBothFixedZero"
            %other_spread=theta_parametric(end-2); w=theta_parametric(end-1);
            %prior = @(s) ((1-w).*exp(-((s-0).^2)./(2.*(sigma_s.^2))))./(sqrt(2*pi).*sigma_s) + (w.*exp(-((s).^2)./(2.*(other_spread.^2))))./(sqrt(2*pi).*other_spread);
            prior = @(s) (1-w).*normpdf(s, 0, sigma_s) + w.* normpdf(s, 0, other_spread);

        case "GaussianLaplaceBothFixedZero"
            %other_spread=theta_parametric(end-2); w=theta_parametric(end-1);
            %prior = @(s) ((1-w).*exp(-((s-0).^2)./(2.*(sigma_s.^2))))./(sqrt(2*pi).*sigma_s) + (w.*exp(-abs(s)./other_spread))./(2.*other_spread);
            prior = @(s) (1-w).*normpdf(s, 0, sigma_s) + w./(2.*other_spread).*exp(-abs(s-0)./other_spread);       

    end
    fun_sum = @(s)  prior(s) .* 1./(sqrt(2.*pi).*sigma_fun(s,sigma0,k).*scale) .* exp(-((x_obs_vals-s).^2)./(2.*((sigma_fun(s,sigma0,k).*scale).^2)));
    fun = @(s) s .* fun_sum(s);
    
    switch numerical_method
        case "midpoint"
            s_hats_PM = midpoint_estimate_integral(fun,a,b,numbins)./midpoint_estimate_integral(fun_sum,a,b,numbins);
        case "trapezoid"
            s_hats_PM = trepazoid_estimate_integral(fun,a,b,numbins)./trepazoid_estimate_integral(fun_sum,a,b,numbins);
    end
    tEnd = toc(tStart);
    
    % Sanity check for computing posterior mean, only used for constant noise
    %s_hats_analytical = (x_obs_vals ./ ((sigma0.*scale).^2)) ./ (1./((sigma0.*scale).^2)+ 1./(sigma_s.^2));
    %s_hats_PM
end