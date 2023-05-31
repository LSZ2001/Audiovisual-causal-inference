function [theta_UA] = complete_thetaua_for_ujointfits(theta_joint, UA_param_keep_idx, UV_has_rescale)
    % theta_joint is in fact [theta_UV, theta_UA_unique];
    % Assume same prior/lapse/sigma_motor for UA and UV. Also same sensory
    % noise type.
    % Only difference is the noise parameters for UV and UA, reflected in
    % theta_UB_unique.
    
    theta_UV = theta_joint(1:(end-length(UA_param_keep_idx)));
    theta_UA_unique = theta_joint((end-length(UA_param_keep_idx)+1):end);
    
    % Remove scale1, scale2, because they are for UV only.
    if(any(UA_param_keep_idx==3)) % there is a k2 -> exp sensory noise model
        rel_scale_idx = [4,5];
    else
        rel_scale_idx = [3,4];
    end
    
    UA_keep_idx = setdiff(1:length(theta_UV), rel_scale_idx);
    theta_UA = theta_UV(UA_keep_idx);
    
    % Add in the rescale for auditory
    if(UV_has_rescale)
        theta_UA(UA_param_keep_idx) = theta_UA_unique; % If UV has rescale, just replace.
    else % If UV does not have rescale, need to add in rescale unique to UA:
        % [sigma0, ks, sigma_s, lapse, (rescale), ...]
        theta_UA = [theta_UA(1:(max(rel_scale_idx))), theta_UA_unique(end), theta_UA((max(rel_scale_idx)+1):end)];
        theta_UA(UA_param_keep_idx) = theta_UA_unique;
    end
    
end