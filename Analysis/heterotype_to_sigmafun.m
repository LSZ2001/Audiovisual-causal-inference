function [sigma_fun] = heterotype_to_sigmafun(hetero_type)
    switch hetero_type
        case "constant"
            sigma_fun = @(s, sigma0, k) sigma0;
        case "linear"
            sigma_fun = @(s, sigma0, k) sigma0 + k(1).*abs(s);
        case "quad"
            sigma_fun = @(s, sigma0, k) sigma0 + k(1).*s.^2;
        case "squarequad"
            sigma_fun = @(s, sigma0, k) sqrt(sigma0.^2 + (k(1).*s).^2);
        case "exp"
            sigma_fun = @(s, sigma0, k) sigma0 + k(1).*(1-exp(-k(2)*abs(s)));
        case "cos"
            sigma_fun = @(s, sigma0, k) sigma0.*(1 + (2.*(k(1).*90./pi).^2) .* (1-cos(s.*pi./90)));
    end
end