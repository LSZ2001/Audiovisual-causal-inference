% This function is to create BADS bounds for Unimodal joint fits of UA and
% UV. It considers which parameters are unique to UA, and concatenates the
% corresponding (unique, unshared) bounds after the bounds of UV data. 

% It assumes identical prior/lapse/sigma_motor for UA and UV, as well as
% same sensory noise type. Only the noise parameters can differ between UA
% and UV.

function [LB, UB, PLB, PUB, UA_param_keep_idx] = merge_UJoint_BADsbounds(LB_UV,UB_UV,PLB_UV,PUB_UV,LB_UA,UB_UA,PLB_UA,PUB_UA, ModelComponents_UA)
    % Merge the BADS bounds of UV and UA to one vector, UV in front.
    if(ModelComponents_UA.SensoryNoise ~= "exp")
        UA_param_keep_idx = [1,2,5]; % sigma0, k, rescale
    else
        UA_param_keep_idx = [1,2,3,6]; % sigma0, k1, k2, rescale
    end
    if(ModelComponents_UA.Rescale~="free") % remove the aud rescale
        UA_param_keep_idx = UA_param_keep_idx(1:(end-1));
    end

    LB = [LB_UV, LB_UA(UA_param_keep_idx)]; % sigma0, ks, rescale
    UB = [UB_UV, UB_UA(UA_param_keep_idx)];
    PLB = [PLB_UV, PLB_UA(UA_param_keep_idx)];
    PUB = [PUB_UV, PUB_UA(UA_param_keep_idx)];
    
%     if(ModelComponents_UA.Rescale~="free") % remove the aud rescale
%         rescale_aud = str2num(ModelComponents_UA.Rescale)
%         LB = [LB, rescale_aud];
%         UB = [UB, rescale_aud];
%         PLB = [PLB, rescale_aud];
%         PUB = [PUB, rescale_aud];
%     end

end