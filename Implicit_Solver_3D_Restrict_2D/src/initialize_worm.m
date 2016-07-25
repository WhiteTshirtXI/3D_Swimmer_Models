function X = initialize_worm(kappa,ds);

  % record the number of points based on the input curvature
  %
  Ns = length(kappa);
    
  % allocate space for position and tangent angle
  %
  X     = zeros(Ns,2);
  th    = zeros(Ns,1);

  % initialize so that s=0 is at the origin and flat
  %
  X(1,1)   = 0;
  X(1,2)   = 0;
  th(1)    = 0;

  % integrate curvature to get tangent angle
  % integrate tangent angle to define position
  %
  for i=2:Ns
    th(i)  = th(i-1)  - ds*kappa(i);
    X(i,1) = X(i-1,1) + ds*cos(th(i-1));
    X(i,2) = X(i-1,2) + ds*sin(th(i-1)); 
  end

  % put the center-of-mass at the origin
  % orient so that the tip-to-tip angle is zero
  %  both of these choices are arbitrary and can be changed externally
  %
  X = X - repmat(mean(X),Ns,1);
  theta = atan2(X(Ns,2)-X(1,2),X(Ns,1)-X(1,1));
  R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
  X = X*R;

  
