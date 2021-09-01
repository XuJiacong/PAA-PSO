function [gbest,gbestval,hist_v,hist_p,fitcount]= PID_CLPSO_func(fhd,Dimension,Particle_Number,Max_Gen,VRmin,VRmax,varargin)
    
    alpha = 0.9;
    kp = 0.08;
    ki = 0.96;
    kd = 0.0;
    
    rand('state',sum(100*clock));
    me = Max_Gen;
    ps = Particle_Number;
    D = Dimension;
    t = 0 : 1/(ps-1) : 1;
    Pc = 0.05 + 0.45.*(exp(10.*t)-1) ./ (exp(10)-1); % or you can change 5 to 10 
    iwt = 0.9 - (1:double(me)) * (0.5/double(me));
    c = 1.49445;
    gapM = 5;
    
    der = zeros(ps, D);
    pro = zeros(ps, D);
    vv_last = zeros(ps, D);
    vv_last2 = zeros(ps, D);
    % Initiate history of gbestvalue
    hist_v = zeros(me, 1);
    % Initiate history of gbest
    hist_p = zeros(me, D);
    
    if length(VRmin)==1
        VRmin=repmat(VRmin,1,D);
        VRmax=repmat(VRmax,1,D);
    end
    
    mv = 0.2 * (VRmax-VRmin);
    Vmin = repmat(-mv, ps, 1);
    Vmax = -Vmin;
    VRmin = repmat(VRmin, ps, 1);
    VRmax = repmat(VRmax, ps, 1);
    
    pos = VRmin + (VRmax-VRmin) .* rand(ps, D);
    e = zeros(ps, 1);
    for i = 1 : ps
        e(i,1) = feval(fhd, pos(i,:)', varargin{:});
    end
    fitcount = ps;
    
    % initialize the velocity of the particles
    vel = Vmin + 2 .* Vmax .* rand(ps, D);
    % initialize the pbest and the pbest's fitness value
    pbest = pos;
    pbestval = e; 
    [gbestval, gbestid] = min(pbestval);
    % initialize the gbest
    gbest = pbest(gbestid, :);
    %gbestrep = repmat(gbest, ps, 1);
    hist_p(1,:) = gbest;
    hist_v(1) = gbestval;

    stay_num = zeros(ps, 1);
    pbest_f = zeros(ps, D);
    
    f_pbest = 1 : ps;
    f_pbest = repmat(f_pbest', 1, D); 

    for k = 1 : ps
        % For each dimension (whose random number is larger than Pc(k)),
        % select the better pbest to learn
        fi1 = ceil(ps*rand(1, D));
        fi2 = ceil(ps*rand(1, D));
        fi = (pbestval(fi1) < pbestval(fi2))' .* fi1 + (pbestval(fi1) >= pbestval(fi2))' .* fi2;
        bi = ceil(rand(1, D) - 1 + Pc(k));
        % If all exemplars of a particle are its own pbest, we will randomly 
        % choose one dimension to learn from another particle’s pbest's 
        % corresponding dimession
        if bi == zeros(1,D)
            rc = randperm(D);
            bi(rc(1)) = 1;
        end
        f_pbest(k, :) = bi.*fi + (1-bi).*f_pbest(k,:);
    end
    

    for i = 2 : me
        for k = 1 : ps
            if stay_num(k) >= gapM
                stay_num(k) = 0;
                f_pbest(k,:) = k .* ones(1, D);
                fi1 = ceil(ps*rand(1, D));
                fi2 = ceil(ps*rand(1, D));
                fi = (pbestval(fi1) < pbestval(fi2))' .* fi1 + (pbestval(fi1) >= pbestval(fi2))' .* fi2;
                bi = ceil(rand(1,D) - 1 + Pc(k));
                if bi == zeros(1,D)
                    rc = randperm(D);
                    bi(rc(1))=1;
                end
                f_pbest(k,:) = bi.*fi + (1-bi).*f_pbest(k,:);
            end

            for dimcnt = 1 : D
                pbest_f(k, dimcnt) = pbest(f_pbest(k, dimcnt), dimcnt);
            end
            
            % Updating Rule
            aa        = c .* rand(1,D) .* (pbest_f(k,:) - pos(k,:));
            vv        = aa;
            vel(k, :) = iwt(i) .* vel(k,:) + aa;
            
            der(k, :) = alpha.*der(k, :) + (vv - 2.*vv_last(k,:) + vv_last2(k,:));
            vv_last2(k,:) = vv_last(k,:);
            
            pro(k, :) = alpha.*pro(k, :) + (vv - vv_last(k,:));
            vv_last(k,:) = vv;
            
            vel(k, :) = (vel(k,:) > mv)   .* mv   + (vel(k,:) <= mv)  .*vel(k,:);
            vel(k, :) = (vel(k,:) < (-mv)).*(-mv) + (vel(k,:) >=(-mv)).*vel(k,:);
            der(k, :) = (der(k,:) > mv)   .* mv   + (der(k,:) <= mv)  .*der(k,:);
            der(k, :) = (der(k,:) < (-mv)).*(-mv) + (der(k,:) >=(-mv)).*der(k,:);
            pro(k, :) = (pro(k,:) > mv)   .* mv   + (pro(k,:) <= mv)  .*pro(k,:);
            pro(k, :) = (pro(k,:) < (-mv)).*(-mv) + (pro(k,:) >=(-mv)).*pro(k,:);
            
            pos(k, :) = pos(k,:) + kp.*pro(k,:) + ki.*vel(k,:) + kd.*der(k,:);
            
            pos(k, :) =   ((pos(k,:)<=VRmax(k,:))&(pos(k,:)>=VRmin(k,:))) .* pos(k,:) ...
                        + (pos(k,:)<VRmin(k,:)) .* (VRmin(k,:)+0.25.*(VRmax(k,:)-VRmin(k,:)).*rand(1,D)) ...
                        + (pos(k,:)>VRmax(k,:)) .* (VRmax(k,:)-0.25.*(VRmax(k,:)-VRmin(k,:)).*rand(1,D));      
            
            e(k,1) = feval(fhd, pos(k,:)', varargin{:});
            fitcount = fitcount + 1;
            
            % Update the pbest
            if e(k,1) < pbestval(k)
                pbest(k,:) = pos(k,:);
                pbestval(k)= e(k);
            else
                stay_num(k) = stay_num(k) + 1;
            end

            % Update the gbest
            if pbestval(k) < gbestval
                gbest = pbest(k, :);
                gbestval = pbestval(k);
                %gbestrep = repmat(gbest, ps, 1);
            end
            hist_p(i,:) = gbest;
            hist_v(i) = gbestval;
        end
    end
end

