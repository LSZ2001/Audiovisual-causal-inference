% Inputs: ModelComponents: struct with fields {
% str         PriorType: "SingleGaussian", "TwoGaussiansOneFixedZero", "TwoGaussiansBothFixedZero"
% str         SensoryNoise: "constant", "linear", "quad", "squarequad", "exp", "cos"
% float/str   Rescale: "free" means to be fitted. Also, can be 1 (no rescale) or 4/3.
% boolean     IsFixedPriorMean: (if true, fit a prior mean parameter called mu; if not, mu is fixed at 0).
% vector      LapseRange: (-45 to 45 for now)}

% theta = [sigma0, k[...], sigma_s, mu?, lapse, rescale?, sigma_center?, w?] is a set of free params values raised by BADS.
%    for hetero_type="constant", still need to provide k (any value makes no difference, but not fitted.
%       the output result for the k column would be 0.
%    mu present only when ~IsFixedPriorMean (hence need ~IsFixedPriorMean for TwoGaussianOneFixedZero)
%    rescale only present when ~FixedRescale && IsFitRescale. 
%    (sigma_center, w) only present for PriorType="TwoGaussiansOneFixedZero".

function [LB, UB, PLB, PUB] = sigmafun_badsbounds_comprehensive(ModelComponents)
    prior_type = ModelComponents.PriorType;
    hetero_type = ModelComponents.SensoryNoise;
    rescale = ModelComponents.Rescale;
    is_fixed_prior_mean = ModelComponents.IsFixedPriorMean;
    is_ML_model = ModelComponents.IsMLModel;
    motor_noise_type =  ModelComponents.MotorNoise;
    if(~isfield(ModelComponents, "NumReliabilityLevels")) % For auditory info, no reliability levels by default.
        num_reliability_levels = 1;
    else
        num_reliability_levels = ModelComponents.NumReliabilityLevels;
    end
    
    if(is_ML_model) % ML model
        LB =  [0.1 0   0   ];
        UB =  [10  3   1   ];
        PLB = [2   0.1 0.01];
        PUB = [5   1   0.2 ];
    elseif(is_fixed_prior_mean || prior_type=="BothGaussiansFixedAtZero") % PM model, additionally fit sigma_s
        LB =  [0.1 0   0.1 0   ];
        UB =  [10  3   30  1   ];
        PLB = [2   0.1 1   0.01];
        PUB = [5   1   10  0.2 ];
    else % PMmufree model, additionally fit mu
        LB =  [0.1 0   0.1 -20 0];
        UB =  [10  3   30  20  1   ];
        PLB = [2   0.1 7   -5  0.01];
        PUB = [5   0.5 10  5   0.2 ];
    end

    k_idx = [2];
    switch hetero_type
        case "constant" % keep 0 as placeholder for k value.
            LB = [LB(1), 0,LB(3:end)];
            UB = [UB(1), 0,UB(3:end)];
            PLB = [PLB(1), 0,PLB(3:end)];
            PUB = [PUB(1), 0,PUB(3:end)];
        case "linear"
        case "quad"
        case "squarequad" % baseline -> hence parameters are above.
        case "exp"
            LB = [LB(1), 0, 0, LB(3:end)];
            UB = [UB(1), 20, 10, UB(3:end)];
            PLB = [PLB(1), 0.001, 0.001, PLB(3:end)];
            PUB = [PUB(1), 10, 3, PUB(3:end)];
            k_idx = [2,3];
        case "cos"
    end
    
        
    if(num_reliability_levels > 1) % params are the multiplication factors for higher reliabilties compared to baseline reliability.
        LB = [LB(1:k_idx(end)), repmat(0.7,1,num_reliability_levels-1), LB((k_idx(end)+1):end)];
        UB = [UB(1:k_idx(end)), repmat(20,1,num_reliability_levels-1), UB((k_idx(end)+1):end)];
        PLB = [PLB(1:k_idx(end)), repmat(1,1,num_reliability_levels-1), PLB((k_idx(end)+1):end)];
        PUB = [PUB(1:k_idx(end)), repmat(5,1,num_reliability_levels-1), PUB((k_idx(end)+1):end)];
    end
    
    if(rescale=="free") % LB/UB/PLB/PUB for rescale param.
        LB = [LB, 1/3];
        UB = [UB, 3];
        PLB = [PLB, 2/3];
        PUB = [PUB, 3/2];
    end
    
    if(prior_type=="TwoGaussiansOneFixedZero") % for (sigma_center, w).
        LB = [LB, 0.1, 0];
        UB = [UB, 30, 1];
        PLB = [PLB, 1, 0.001];
        PUB = [PUB, 10, 0.999];
    elseif(prior_type=="TwoGaussiansBothFixedZero") % for (delta_sigma, w).
        % Note: here the 2 Gaussian prior mixture components are
        % interchangeable. To avoid confusion, we have them  sigma_1 =
        % simga_s, and sigma2 = sigma_s + delta_sigma.
        LB = [LB, 0, 0];
        UB = [UB, 20, 1];
        PLB = [PLB, 1, 0.001];
        PUB = [PUB, 10, 0.999];
    elseif(prior_type=="GaussianLaplaceBothFixedZero")
        LB = [LB, 0.1, 0];
        UB = [UB, 20, 1];
        PLB = [PLB, 1, 0.001];
        PUB = [PUB, 10, 0.999];
        
    end
    
    if(motor_noise_type=="Gaussian") % decision noise is centered at s_PM = F(x_obs), with std to be fitted.
        LB = [LB, 0.25]; % LB of motor_noise must >> dx of x_grid used in NLL computation code. 
        UB = [UB, 2];
        PLB = [PLB, 0.5];
        PUB = [PUB, 1.5];
    end
        
        
end