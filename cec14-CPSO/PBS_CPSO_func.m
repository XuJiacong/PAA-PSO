function [gbest,gbestval,hist_v,hist_p,fitcount]= PBS_CPSO_func(fhd,Dimension,Particle_Number,iter_max, VRmin,VRmax,varargin)
    % CPSO-H6
    kd = 0.03;
    alpha = 0.9;
    ps = Particle_Number;% Particle Number
    D = Dimension; % Dimension of the Target Problem
    cc = [1.49 1.49]; % Acceleration Constants
    groupnum = 6; % Number of Groups
    %me = ceil(Max_FES / ps / 6); % MAX Generation
    me = iter_max;
    iwt = 0.9 - (1:me) .* (0.5 ./ me);% Inertia Weight
    e1 = zeros(ps,1);
    e = zeros(ps,1);
    vel = zeros(ps,D);
    vel1 = zeros(ps,D);
    aa = zeros(ps,D);
    aa1 = zeros(ps,D);
    der = zeros(ps,D);
    last_gbestrep = zeros(ps,D);
    % Initiate history of gbestvalue
    hist_v = zeros(me, 1);
    % Initiate history of gbest
    hist_p = zeros(me, D);

%% Boundaries

    if length(VRmin) == 1
        VRmin = repmat(VRmin,1,D);
        VRmax = repmat(VRmax,1,D);
    end
    mv = 0.2 * (VRmax - VRmin);
    mv1 = mv;
    VRmin = repmat(VRmin,ps,1);
    VRmax = repmat(VRmax,ps,1);

%% Initialize Every Particle

    pos = VRmin + (VRmax - VRmin) .* rand(ps,D);
    pos1 = VRmin + (VRmax - VRmin) .* rand(ps,D);

%% First Loop

    for j = 1:ps
        e(j,1) = feval(fhd,pos(j,:)',varargin{:});
        e1(j,1) = feval(fhd,pos1(j,:)',varargin{:});
    end
    fitcount = ps;
    %initialize the pbest and the pbest's fitness value
    pbest = pos;
    pbest1 = pos1;
    pbestval = e;
    pbestval1 = e1;
    
    %initialize the gbest and the gbest's fitness value
    [gbestval,gbestid] = min(pbestval);
    gbest = pbest(gbestid,:);
    hist_v(1) = gbestval;
    hist_p(1,:) = gbest;
    gbestrep = repmat(gbest,ps,1);
    [gbestval1,gbestid1] = min(pbestval1);
    gbest1 = pbest1(gbestid1,:);
    gbestrep1 = repmat(gbest1,ps,1);
    gbestid1 = repmat(gbestid1,1,groupnum);
    pbestval1 = repmat(pbestval1,1,groupnum);
    e1 = repmat(e1,1,groupnum);

%% Assign Group Members

    rc = randperm(D);
    tmp = round(D/groupnum);
    particle_num(1,:) = rc(1:tmp);
    for j = 1:groupnum-1
        particle_num(j,:) = rc((j*tmp-tmp+1):(j*tmp));
    end
    particle_num(groupnum,:) = rc((D-tmp+1):D);
    group_d = particle_num;

%% Main Loop

    for i = 2:me
        for j = 1:groupnum
            len_grpd = length(group_d(j,:));
            for k = 1:ps
                pos1_temp = gbest1;
                pos1_temp(group_d(j,:)) = pos1(k,group_d(j,:));

                e1(k,j) = feval(fhd,pos1_temp',varargin{:});
                fitcount = fitcount + 1;
                
                if e1(k,j) < pbestval1(k,j)
                    pbest1(k,group_d(j,:)) = pos1(k,group_d(j,:));
                    pbestval1(k,j) = e1(k,j);
                end

                if pbestval1(k,j) < gbestval1
                    gbest1(group_d(j,:)) = pbest1(k,group_d(j,:));
                    gbestval1 = pbestval1(k,j);
                    gbestrep1 = repmat(gbest1,ps,1);
                    gbestid1(j) = k;
                end

                % Update Acceleration Vector
                aa1(k,group_d(j,:)) = cc(1) .* rand(1, len_grpd) .* (pbest1(k,group_d(j,:)) - ... 
                                      pos1(k,group_d(j,:))) + cc(2) .* rand(1, len_grpd) .* ...
                                      (gbestrep1(k,group_d(j,:)) - pos1(k,group_d(j,:)));

                % Update Velocity Vector
                vel1(k,group_d(j,:)) = iwt(i) .* vel1(k,group_d(j,:)) + aa1(k,group_d(j,:));
                % Clamping
                vel1(k,group_d(j,:)) = (vel1(k,group_d(j,:)) > mv1(group_d(j,:))) .* mv1(group_d(j,:)) + ...
                                       (vel1(k,group_d(j,:)) <= mv1(group_d(j,:))) .* vel1(k,group_d(j,:));
                vel1(k,group_d(j,:)) = (vel1(k,group_d(j,:)) < (-mv1(group_d(j,:)))) .* (-mv1(group_d(j,:))) + ...
                                       (vel1(k,group_d(j,:)) >= (-mv1(group_d(j,:)))) .* vel1(k,group_d(j,:));

                % Update Position
                pos1(k,group_d(j,:)) = pos1(k,group_d(j,:)) + vel1(k,group_d(j,:));
                
                pos1(k,group_d(j,:)) = (pos1(k,group_d(j,:)) > VRmax(k,group_d(j,:))) .* (VRmax(k,group_d(j,:))-0.25*rand(1,len_grpd).*(VRmax(k,group_d(j,:))-VRmin(k,group_d(j,:)))) + ...
                                       (pos1(k,group_d(j,:)) <= VRmax(k,group_d(j,:))) .* pos1(k,group_d(j,:));
                pos1(k,group_d(j,:)) = (pos1(k,group_d(j,:)) < VRmin(k,group_d(j,:))) .* (VRmin(k,group_d(j,:))+0.25*rand(1,len_grpd).*(VRmax(k,group_d(j,:))-VRmin(k,group_d(j,:)))) + ...
                                       (pos1(k,group_d(j,:)) >= VRmin(k,group_d(j,:))) .* pos1(k,group_d(j,:));

            end
        end
        
        % Replace a random agent with the best one in 1 
        rc = randperm(ps);
        k = rc(1);
        if k == gbestid
            k = rc(2);
        end
        pos(k,:) = gbest1;

        for k = 1:ps
            e(k,1) = feval(fhd,pos(k,:)',varargin{:});
            fitcount = fitcount + 1;
            % Update the pbest
            if  e(k) < pbestval(k)
                pbest(k,:) = pos(k,:);
                pbestval(k) = e(k);
            end
            % Update the gbest
            if pbestval(k) < gbestval
                gbest = pbest(k,:);
                gbestid = k;
                gbestval = pbestval(k);
                gbestrep = repmat(gbest,ps,1);
            end

            % Update Acceleration Vector
            aa(k,:) = cc(1) .* rand(1,D) .* (pbest(k,:) - pos(k,:)) + ...
                      cc(2) .* rand(1,D) .* (gbestrep(k,:) - pos(k,:));
            % Update der vector
            der(k,:) = alpha*der(k,:) + (gbestrep(k,:) - last_gbestrep(k,:));
            last_gbestrep(k,:) = gbestrep(k,:);
            % Update Velocity Vector
            vel(k,:) = iwt(i) .* vel(k,:) + aa(k,:);
            
            vel(k,:) = (vel(k,:) > mv) .* mv + (vel(k,:) <= mv) .* vel(k,:);
            vel(k,:) = (vel(k,:) < (-mv)) .* (-mv) + (vel(k,:) >= (-mv)) .* vel(k,:);

            % Update Position of every Particl
            pos(k,:) = pos(k,:) + vel(k,:) + kd*der(k,:);

            pos(k,:) = (pos(k,:) > VRmax(k,:)) .* (VRmax(k,:)-0.25.*(VRmax(k,:)-VRmin(k,:)).*rand(1,D)) + ...
                       (pos(k,:) <= VRmax(k,:)) .* pos(k,:);
            pos(k,:) = (pos(k,:) < VRmin(k,:)) .* (VRmin(k,:)+0.25.*(VRmax(k,:)-VRmin(k,:)).*rand(1,D)) + ...
                       (pos(k,:) >= VRmin(k,:)) .* pos(k,:);

        end

        for j = 1:groupnum
            rc = randperm(ps);
            k = rc(1);
            if k == gbestid1(j)
                k = rc(2);
            end
            pos1(k,group_d(j,:)) = gbest(group_d(j,:));
        end
    hist_v(i) = gbestval;
    hist_p(i,:) = gbest;
    end
end