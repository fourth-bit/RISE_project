function out = wf_inv_exp(trigger_amount)
    trigger_amount(trigger_amount == 1) = 0;
    out = exp(-4 * trigger_amount) - exp(-4);
end

