function [out] = wf_none_orig(trigger_amount)
% The default all-or-nothing behavior that is in the original paper
out = trigger_amount == 1;
end

