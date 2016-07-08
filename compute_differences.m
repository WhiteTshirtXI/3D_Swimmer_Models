
% output locations
%
datadir    = './diffdata';  
% time stepping info
%
t0    = 0.1;
dtout = 0.1;
Tend  = 2.0;

% record the number of outputs of position
%
k=1;
for t = t0:dtout:Tend
  filename1 = sprintf('./data/imworm_n064_t%f.mat',t);
  data1 = load(filename1);
  filename2 = sprintf('./data/exworm_n064_t%f.mat',t);
  data2 = load(filename2);
  foutn1 = sprintf('%s/%s_t%f.mat',datadir,'fluidVelocity',t);
  foutn2 = sprintf('%s/%s_t%f.mat',datadir,'wormVelocity',t);
  foutn3 = sprintf('%s/%s_t%f.mat',datadir,'position',t);
  foutn4 = sprintf('%s/%s_t%f.mat',datadir,'stress',t);
  foutn5 = sprintf('%s/%s_t%f.mat',datadir,'speed',t);
  speed1 = get_speed(data1.XTworm, data1.Nt, data1.Tend, data1.Tper);
  speed2 = get_speed(data2.XTworm, data2.Nt, data2.Tend, data2.Tper);
  fluidForce = absDifference(data1.U, data2.U);
  wormForce = absDifference(data1.Uw, data2.Uw);
  wormPosition = absDifference(data1.XTworm, data2.XTworm);
  stressForce = absDifference(data1.Shat, data2.Shat);
  speed = absDifference(speed1, speed2);
  save(foutn1, 'fluidForce');
  save(foutn2, 'wormForce');
  save(foutn3, 'wormPosition');
  save(foutn4, 'stressForce');
  save(foutn5, 'speed');
end
