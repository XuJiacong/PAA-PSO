function [gbest, gbestval, hist_v, hist_p, fitcount]= PID_PSO_funcM(fhd, Dimension, Particle_Number, Max_iters, VRmin, VRmax, varargin)
    gamma = (0.99 - (1:double(Max_iters)) .* (0.98./double(Max_iters)));
    alpha = 0.9;
    kp = 0.05;
    %ki = 1.00;
    %kd = 0.0;
    % Random initialization
    rand('state', sum(100*clock));
    % Number of particles
    Num_Ptcs = Particle_Number;
    % Number of dimensions
    D = Dimension;
    % Acceleration constants
    C = [2 2];
    % Inertia weight decreases linearly from 0.9 to 0.4
    Iwt = (0.9 - (1:double(Max_iters)) .* (0.5./double(Max_iters)));
    
    % Initiate history of gbestvalue
    hist_v = zeros(Max_iters, 1);
    % Initiate history of gbest
    hist_p = zeros(Max_iters, D);
    
    % The bounds of variables
    if length(VRmin) == 1
        VRmin = repmat(VRmin, 1, D);
        VRmax = repmat(VRmax, 1, D);
    end
    % The bounds of velocities
    mv = 0.2 * (VRmax - VRmin);
    Vmin = repmat(-mv, Num_Ptcs, 1);
    Vmax = -Vmin;
    % Extended to all particles
    VRmin = repmat(VRmin, Num_Ptcs, 1);
    VRmax = repmat(VRmax, Num_Ptcs, 1);
    % Randomly initiate positions of particles 
    Pos = VRmin + (VRmax-VRmin) .* rand(Num_Ptcs, D);
    % Evaluate the error value for each initiated position
    e = feval(fhd, Pos', varargin{:}); % 1 x Num_Ptcs
    % Fitness evaluations add number of particles each time
    fitcount = Num_Ptcs;
    % Initialize the velocity of the particles
    vel = Vmin + 2 .* Vmax .* rand(Num_Ptcs, D);
    % Initialize the pbest and the pbest's fitness value
    pbest = Pos;
    pbestval = e;
    % Initialize the gbest and the gbest's fitness value
    [gbestval,gbestid] = min(pbestval);
    hist_v(1) = gbestval;
    gbest = pbest(gbestid, :);
    hist_p(1,:) = gbest;
    gbestrep = repmat(gbest, Num_Ptcs, 1);
    gbestrep_last = zeros(Num_Ptcs, D);
    pbest_last = zeros(Num_Ptcs, D);
%     gbestrep_last2 = zeros(Num_Ptcs, D);
%     pbest_last2 = zeros(Num_Ptcs, D);
%     der = zeros(Num_Ptcs, D);
    pro = zeros(Num_Ptcs, D);
    
    for i = 2 : Max_iters
        
        r1 = rand(Num_Ptcs, D);
        r2 = rand(Num_Ptcs, D);
        vv = C(1) .* r1 .* (pbest-Pos) + C(2) .* r2 .* (gbestrep-Pos);
        
        % Update D component
%         der = alpha.*der + (gamma(i)) .* (pbest-2.*pbest_last+pbest_last2) + (1-gamma(i)) .* (gbestrep-2.*gbestrep_last+gbestrep_last2);
%         pbest_last2 = pbest_last;
%         gbestrep_last2 = gbestrep_last;
        
        % Update P component
        pro = alpha.*pro + (gamma(i)) .* (pbest-pbest_last) + (1-gamma(i)) .* (gbestrep-gbestrep_last);
        pbest_last = pbest;
        gbestrep_last = gbestrep;
        
        % Update velocity (I component)
        vel = Iwt(i) .* vel + vv;
        
        % Velocity clamping
        vel = (vel>Vmax) .* Vmax + (vel<=Vmax) .* vel;
        vel = (vel<Vmin) .* Vmin + (vel>=Vmin) .* vel;
              
        % Update positions
        Pos = Pos + kp.*pro + vel;% + kd.*der;
        
        % Return the particles outside the bounds
        Pos = ((Pos>=VRmin)&(Pos<=VRmax)) .*Pos...
              + (Pos<VRmin) .*(VRmin+0.25.*(VRmax-VRmin).*rand(Num_Ptcs,D))...
              + (Pos>VRmax) .*(VRmax-0.25.*(VRmax-VRmin).*rand(Num_Ptcs,D));
          
        % Fitness evaluation
        e = feval(fhd, Pos', varargin{:});
        fitcount = fitcount + Num_Ptcs;
        
        % Update the pbest
        tmp = pbestval < e;
        temp = repmat(tmp', 1, D);
        pbest = temp.*pbest + (1-temp).*Pos;
        pbestval = tmp.*pbestval + (1-tmp).*e;
        
        % Update the gbest
        [gbestval,tmp] = min(pbestval);
        hist_v(i) = gbestval;
        gbest = pbest(tmp, :);
        hist_p(i,:) = gbest;
        gbestrep = repmat(gbest, Num_Ptcs, 1);
    end
end