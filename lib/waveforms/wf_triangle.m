function [out] = wf_triangle(trigger_amount)
    out = clip(-4 * abs(trigger_amount - 0.25) + 1, 0, 1);
end

