function [ave] = getPeriodAverage(XTworm, Tperstep)
endStep = size(XTworm,3);
[MAV,U,Xcm] = get_speed(XTworm(:,:,endStep-Tperstep:endStep));
ave = sum(MAV)/size(MAV,1);
end

