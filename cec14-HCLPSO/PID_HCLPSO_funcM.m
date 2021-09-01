function [gbest,gbestval,hist_v,hist_p,FES]= PID_HCLPSO_funcM(fhd,dimension,num_particle,max_iters,VRmin,VRmax,varargin) 

kd = 0.05;
alpha = 0.9;
rand('state',sum(100*clock));
ps = num_particle;
me = max_iters;
D = dimension;
gamma = (0.9 - (1:double(me)) .* (0.8./double(me)));
% Initiate history of gbestvalue
hist_v = zeros(me, 1);
% Initiate history of gbest
hist_p = zeros(me, D);
der = zeros(ps,D);

last_gbest = zeros(ps,D);
last_pbest = zeros(ps,D);

num_g1=15;
num_g2=ps-num_g1;

j=0:(1/(ps-1)):1; % Learning Probability Curve Pc
j=j*10;
Pc=ones(D,1)*(0.0+((0.25).*(exp(j)-exp(j(1)))./(exp(j(ps))-exp(j(1))))); 

Weight =0.99-(1:me)*0.79/me; % Inertia Weight

K=3-(1:me)*1.5/me;% Acceleration Coefficients
c1=2.5-(1:me)*2/me;
c2=0.5+(1:me)*2/me;

% Initialization
interval = VRmax-VRmin;
Vmax=interval*0.2; 
Vmin=-Vmax;
pos = VRmin+ (VRmax-VRmin).*rand(ps,D);    %  position
vel = Vmin + (Vmax-Vmin).*rand(ps,D); % velocity

e = feval(fhd, pos', varargin{:});
fitcount=ps;

[gbestval,g_index]=min(e);
gbest=pos(g_index,:);
hist_p(1,:) = gbest;
hist_v(1) = gbestval;

pbest=pos; 
pbestval=e;

obj_func_slope=zeros(ps,1);
fri_best=(1:ps)'*ones(1,D);
fri_best_pos = zeros(ps,D);
        
for i=1:num_g1 % Updateding examplers for; group 1
    fri_best(i,:)=i*ones(1,D);
    friend1=ceil(num_g1*rand(1,D));
    friend2=ceil(num_g1*rand(1,D));
    friend=(pbestval(friend1)<pbestval(friend2)).*friend1+(pbestval(friend1)>=pbestval(friend2)).*friend2;
    toss=ceil(rand(1,D)-Pc(:,i)');
    if toss==ones(1,D)
       temp_index=randperm(D);
       toss(1,temp_index(1))=0;
       clear temp_index;
    end
    fri_best(i,:)=(1-toss).*friend+toss.*fri_best(i,:);
    for d=1:D
        fri_best_pos(i,d)=pbest(fri_best(i,d),d);
    end
end

for i=num_g1+1:ps % Updateding examplers for; group 1
    fri_best(i,:)=i*ones(1,D);
    friend1=ceil(ps*rand(1,D));
    friend2=ceil(ps*rand(1,D));
    friend=(pbestval(friend1)<pbestval(friend2)).*friend1+(pbestval(friend1)>=pbestval(friend2)).*friend2;
    toss=ceil(rand(1,D)-Pc(:,i)');
    if toss==ones(1,D)
       temp_index=randperm(D);
       toss(1,temp_index(1))=0;
       clear temp_index;
    end
    fri_best(i,:)=(1-toss).*friend+toss.*fri_best(i,:);
    for d=1:D
        fri_best_pos(i,d)=pbest(fri_best(i,d),d);
    end
end


for k = 1:me-1
        gbestrep=repmat(gbest,ps,1);
        % group 1
         
        delta_g1=(K(k).*rand(num_g1,D).*(fri_best_pos(1:num_g1,:)-pos(1:num_g1,:)));
        vel_g1=Weight(k)*vel(1:num_g1,:)+delta_g1;
        der(1:num_g1,:) = alpha*der(1:num_g1,:) + gamma(k)*(pbest(1:num_g1,:)-last_pbest(1:num_g1,:))+ (1-gamma(k))*(gbestrep(1:num_g1,:)-last_gbest(1:num_g1,:));
        last_pbest(1:num_g1,:) = pbest(1:num_g1,:);
        last_gbest(1:num_g1,:) = gbestrep(1:num_g1,:);
        vel_g1=((vel_g1<Vmin).*Vmin)+((vel_g1>Vmax).*Vmax)+(((vel_g1<Vmax)&(vel_g1>Vmin)).*vel_g1);
        pos_g1=pos(1:num_g1,:)+vel_g1+kd*der(1:num_g1,:);

       % group 2 (Only apply in group 2)
        gbest_pos_temp=repmat(gbest,num_g2,1);
        delta_g2=(c1(k).*rand(num_g2,D).*(fri_best_pos(num_g1+1:end,:)-pos(num_g1+1:end,:)))+(c2(k).*rand(num_g2,D).*(gbest_pos_temp-pos(num_g1+1:end,:)));
        vel_g2=Weight(k)*vel(num_g1+1:end,:)+delta_g2;
        vel_g2=((vel_g2<Vmin).*Vmin)+((vel_g2>Vmax).*Vmax)+(((vel_g2<Vmax)&(vel_g2>Vmin)).*vel_g2);
        der(num_g1+1:end,:) = alpha*der(num_g1+1:end,:) + gamma(k)*(pbest(num_g1+1:end,:)-last_pbest(num_g1+1:end,:))+(1-gamma(k))*(gbestrep(num_g1+1:end,:)-last_gbest(num_g1+1:end,:));
        %der(num_g1+1:end,:) = alpha*der(num_g1+1:end,:) + (1-gamma(k))*(gbestrep(num_g1+1:end,:)-last_gbest(num_g1+1:end,:));
        last_pbest(num_g1+1:end,:) = pbest(num_g1+1:end,:);
        last_gbest(num_g1+1:end,:) = gbestrep(num_g1+1:end,:);
        pos_g2=pos(num_g1+1:end,:)+vel_g2+kd*der(num_g1+1:end,:);

         % whole group
    
        pos=[pos_g1;pos_g2];
        vel=[vel_g1;vel_g2];
   
        % Evaluate fitness
         for i=1:ps   
            pos(i, :) =   ((pos(i,:)<=VRmax)&(pos(i,:)>=VRmin)) .* pos(i,:) ...
                        + (pos(i,:)<VRmin) .* (VRmin+0.25.*(VRmax-VRmin).*rand(1,D)) ...
                        + (pos(i,:)>VRmax) .* (VRmax-0.25.*(VRmax-VRmin).*rand(1,D));
            e(i)=feval(fhd, pos(i,:)', varargin{:}); 
            fitcount=fitcount+1;

            if  e(i)<pbestval(i) % update pbest value and position
                pbest(i,:)=pos(i,:);   
                pbestval(i)=e(i);
                obj_func_slope(i)=0;
            else
                obj_func_slope(i)=obj_func_slope(i)+1;
            end
            
            if  pbestval(i)<gbestval % update gbest value and postion
                gbest=pbest(i,:); 
                gbestval=pbestval(i);
            end   
        end 
    
        for i=1:num_g1 % updateding exampler for group 1
            if obj_func_slope(i)>5
                fri_best(i,:)=i*ones(1,D); % for its own pbest
                friend1=ceil(num_g1*rand(1,D));
                friend2=ceil(num_g1*rand(1,D));
                friend=(pbestval(friend1)<pbestval(friend2)).*friend1+(pbestval(friend1)>=pbestval(friend2)).*friend2;
                toss=ceil(rand(1,D)-Pc(:,i)');
            
                if toss==ones(1,D)
                    temp_index=randperm(D);
                    toss(1,temp_index(1))=0;
                    clear temp_index;
                end
            
                fri_best(i,:)=(1-toss).*friend+toss.*fri_best(i,:);
                for d=1:D
                    fri_best_pos(i,d)=pbest(fri_best(i,d),d);
                end
                obj_func_slope(i)=0;
            end
        end % updating exampler for group 1
    
        for i=num_g1+1:ps % updating exampler for group 2
            if obj_func_slope(i)>5
                fri_best(i,:)=i*ones(1,D);
                friend1=ceil(ps*rand(1,D));
                friend2=ceil(ps*rand(1,D));
                friend=(pbestval(friend1)<pbestval(friend2)).*friend1+(pbestval(friend1)>=pbestval(friend2)).*friend2;
                toss=ceil(rand(1,D)-Pc(:,i)');
            
                if toss==ones(1,D)
                    temp_index=randperm(D);
                    toss(1,temp_index(1))=0;
                    clear temp_index;
                end
            
                fri_best(i,:)=(1-toss).*friend+toss.*fri_best(i,:);
                for d=1:D
                fri_best_pos(i,d)=pbest(fri_best(i,d),d);
                end
                obj_func_slope(i)=0;
            end
        end % updating exampler for group 2
        
hist_p(k+1,:) = gbest;
hist_v(k+1) = gbestval; 
end

% fprintf('Output\n');
FES=fitcount;

end
