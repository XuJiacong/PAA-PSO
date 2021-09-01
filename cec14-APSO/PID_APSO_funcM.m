function [gbest,gbestval,hist_v,hist_p,fitcount]= PID_APSO_funcM(fhd,Dimension,Particle_Number,Max_Gen,VRmin,VRmax,varargin)
gamma = (0.9 - (1:double(Max_Gen)) .* (0.8./double(Max_Gen)));
c1 = 2; 
c2 = 2;  
alpha = 0.9;
kp = 0.05;
ki = 1.00;
kd = 0.00;
ps = Particle_Number;
D = Dimension; 
pos = VRmin + (VRmax - VRmin).*rand(ps,D);  
Vmax = 0.2 * (VRmax-VRmin); 
Vmin = -Vmax;
vel = Vmin + 2 .* Vmax .* rand(ps, D);

% Initiate history of gbestvalue
hist_v = zeros(Max_Gen, 1);
% Initiate history of gbest
hist_p = zeros(Max_Gen, D);

% Elitist Learning Rate
sigmax=1; 
sigmin=0.1; 
sig=1;  

e = feval(fhd, pos', varargin{:});
fitcount = ps;
[bestfitness, bestindex]=min(e); 
gbest=pos(bestindex,:);   
pbest=pos;    
fitnesspbest=e; 
gbestval=bestfitness; 
hist_p(1,:) = gbest;
hist_v(1) = gbestval;
der = zeros(ps, D);
pro = zeros(ps, D);
gbest_last = zeros(ps, D);
pbest_last = zeros(ps, D);
gbest_last2 = zeros(ps, D);
pbest_last2 = zeros(ps, D);

for i=2:Max_Gen 
    
    factor=calfactor(pos,ps,gbest); 
    if i==2 
        ind_1=1;
    else
        ind_1=ind; 
    end 
    ind=fuzzyclassification(factor,ind_1);
    
    if ind == 1 
        c1=c1+unifrnd(0.05,0.1); 
        c2=c2-unifrnd(0.05,0.1); 
    elseif ind == 2 
        c1=c1+0.5*unifrnd(0.05,0.1); 
        c2=c2-0.5*unifrnd(0.05,0.1); 
    elseif ind == 3 
        c1=c1+0.5*unifrnd(0.05,0.1); 
        c2=c2+0.5*unifrnd(0.05,0.1); 
        p = gbest; 
        d = unidrnd(D); 
        p(d)=p(d)+(VRmax-VRmin)*normrnd(0,sig^2); 
        p(p(:)>VRmax)=VRmax; 
        p(p(:)<VRmin)=VRmin; 
        p_e = feval(fhd, p', varargin{:});
        if p_e < gbestval 
            gbest = p; 
        else 
            [~,bb] = max(e); 
            pos(bb,:) = p; 
        end 
    else 
        c1=c1-unifrnd(0.05,0.1); 
        c2=c2+unifrnd(0.05,0.1); 
    end 
    
    w=1/(1+1.5*exp(-2.6*factor)); 
    
    % Clamping c1 and c2
    if c1<1.5 
        c1=1.5; 
    elseif c1>2.5 
        c1=2.5; 
    end 
    if c2<1.5 
        c2=1.5; 
    elseif c2>2.5 
        c2=2.5; 
    end 
    crange=c1+c2;
    if crange > 4
        c1=(c1/crange)*4; 
        c2=(c2/crange)*4; 
    end
    
    sig=sigmax-(sigmax-sigmin)*((i-1)/(Max_Gen-1));
    
    vel = w*vel + c1.*rand(ps,D).*(pbest - pos) + c2.*rand(ps,D).*(repmat(gbest,ps,1) - pos);
    % Update D component
    gbestrep = repmat(gbest, ps, 1);
    der = alpha.*der + (gamma(i)) .* (pbest-2.*pbest_last+pbest_last2) + (1-gamma(i)) .* (gbestrep-2.*gbest_last+gbest_last2);
    pbest_last2 = pbest_last;
    gbest_last2 = gbest_last;
        
    % Update P component
    pro = alpha.*pro + (gamma(i)) .* (pbest-pbest_last) + (1-gamma(i)) .* (gbestrep-gbest_last);
    pbest_last = pbest;
    gbest_last = gbestrep;
    
    vel=(vel>Vmax).*Vmax+(vel<=Vmax).*vel;
    vel=(vel<Vmin).*Vmin+(vel>=Vmin).*vel;
    
    pos = pos + kp.*pro + ki.*vel + kd.*der;
    
    pos=((pos>=VRmin)&(pos<=VRmax)).*pos...
            +(pos<VRmin).*(VRmin+0.25.*(VRmax-VRmin).*rand(ps,D))+(pos>VRmax).*(VRmax-0.25.*(VRmax-VRmin).*rand(ps,D));
     
    e = feval(fhd, pos', varargin{:});
    fitcount = fitcount + ps;
    
    for j=1:ps 
         
        if e(j) < fitnesspbest(j) 
            pbest(j,:) = pos(j,:); 
            fitnesspbest(j) = e(j); 
        end 
         

        if e(j) < gbestval 
            gbest = pos(j,:); 
            gbestval = e(j); 
        end 
    end 
    hist_p(i,:) = gbest;
    hist_v(i) = gbestval;    
end 

end