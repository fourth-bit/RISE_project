function out = wf_exp(trigger_amount)
%WF_EXP Summary of this function goes here
%   Detailed explanation goes here
out = exp(8 * (trigger_amount - 1)) - exp(-8);
end

