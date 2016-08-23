function [absDiff] = worm_period_difference(XTworm, Tperstep)
endStep = size(XTworm,3);
absDiff = abs(getPeriodAverage(XTworm, Tperstep) - getPeriodAverage(XTworm(:,:,1:endStep - Tperstep), Tperstep));
end

