function [gbest,gbestval,hist_v,hist_p,FES]= PBS_EPSO_func(fhd,Dimension,num_particle,max_iteration,VRmin,VRmax,varargin)
kd = 0.03;
alpha = 0.9;
% Description 
% CLPSO/Sa-PSO(PSO,FDR,HPSO,LIPS,CLPSO-gbest)

% CLPSO: c1=3~1.5, w1=0.9~0.4
% Sa-PSO: PSO/FDR/HPSO/LIPS/CLPSO-gbest, w=0.9~0.2
% PSO: c2_1=2.5~0.5, c2_2=0.5~2.5
% HPSO: c3_1=2.5~0.5, c3_2=0.5~2.5
% LIPS: nsize=3
% CLPSO with gbest:c4_1=2.5~0.5, c4_2=0.5~2.5 & Pc=0.5
rand('state',sum(100*clock));
me = max_iteration;
D = Dimension;

% Initiate history of gbestvalue
hist_v = zeros(me, 1);
% Initiate history of gbest
hist_p = zeros(me, D);

ps = num_particle;
num_g1 = 8;
num_g2 = num_particle - num_g1;
fri_best_pos = zeros(ps,D);
% Initial position
pos = VRmin + (VRmax - VRmin).*rand(ps,D); 
% Initial velocity
Vmax = 0.2 * (VRmax-VRmin); 
Vmin = -Vmax;
vel = Vmin + 2 .* Vmax .* rand(ps, D);
delta = zeros(ps,D);
Pnd = zeros(ps,D);
der = zeros(ps,D);
last_gbest = zeros(ps,D);

e = feval(fhd, pos', varargin{:});
fitcount = ps;

% Initial global best position
[gbestval,gbest_index] = min(e);
gbest = pos(gbest_index,:); 
hist_p(1,:) = gbest;
hist_v(1) = gbestval;
% Initial pbest position
pbest = pos;
pbest_val = e;

% Method: CLPSO
c1=3-(1:me)*1.5/me;
w1=0.9-(1:me)*(0.5/me);

% Sa-PSO (PSO/FDR/HPSO/LIPS/CLPSO-gbest)
w2=0.9-(1:me)*(0.7/me);

% Method 1: PSO
c2_1=2.5-(1:me)*2/me;
c2_2=0.5+(1:me)*2/me;

% Method 2: FDR_PSO
fii=[1 1 2];

% Method 3: HPSO_TVAC
c3_1=2.5-(1:me)*2/me;
c3_2=0.5+(1:me)*2/me;
re_init_vel=VRmax-(1:me)*(0.9*VRmax)/me;

% Method 4: LIPS
nsize=3;

% Method 5: CLPSO
c4_1=2.5-(1:me)*2/me;
c4_2=0.5+(1:me)*2/me;
obj_func_slope=zeros(ps,1);
fri_best=(1:ps)'*ones(1,D);
j=0:(1/(ps-1)):1; 
j=j*10;
Pc=ones(D,1)*(0.0+((0.5).*(exp(j)-exp(j(1)))./(exp(j(ps))-exp(j(1))))); 

for i=1:num_g1
    fri_best(i,:)=i*ones(1,D);
    friend1=ceil(num_g1*rand(1,D));
    friend2=ceil(num_g1*rand(1,D));
    friend=(pbest_val(friend1)<pbest_val(friend2)).*friend1+(pbest_val(friend1)>=pbest_val(friend2)).*friend2;
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

for i=num_g1+1:num_g1+num_g2
    fri_best(i,:)=i*ones(1,D);
    friend1=ceil(ps*rand(1,D));
    friend2=ceil(ps*rand(1,D));
    friend=(pbest_val(friend1)<pbest_val(friend2)).*friend1+(pbest_val(friend1)>=pbest_val(friend2)).*friend2;
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

num_strategy = 5;
LP=50; % Learning Period

for k = 1:me-1
        gbestrep=repmat(gbest,ps,1);
        
if k<=(0.9*me)
       
      % Group 1: CLPSO
      for i=1:num_g1 
            if obj_func_slope(i)>5
                       fri_best(i,:)=i*ones(1,D);
                       friend1=(ceil(num_g1*rand(1,D)));
                       friend2=(ceil(num_g1*rand(1,D)));
                       friend=(pbest_val(friend1)<pbest_val(friend2)).*friend1+(pbest_val(friend1)>=pbest_val(friend2)).*friend2;
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
            end % if obj_func_slope(i)>5    
             
            delta(i,:)=(c1(1).*rand(1,D).*(fri_best_pos(i,:)-pos(i,:)));
            vel(i,:)=w1(k)*vel(i,:)+delta(i,:);
            vel(i,:)=((vel(i,:)<Vmin).*Vmin)+((vel(i,:)>Vmax).*Vmax)+(((vel(i,:)<Vmax)&(vel(i,:)>Vmin)).*vel(i,:));
            %der(i,:)= alpha*der(i,:) + (gbestrep(i,:)-last_gbest(i,:));
            %last_gbest(i,:) = gbestrep(i,:);
            pos(i,:)=pos(i,:)+vel(i,:);%+kd*der(i,:); 
            
            pos(i, :) =   ((pos(i,:)<=VRmax)&(pos(i,:)>=VRmin)) .* pos(i,:) ...
                        + (pos(i,:)<VRmin) .* (VRmin+0.25.*(VRmax-VRmin).*rand(1,D)) ...
                        + (pos(i,:)>VRmax) .* (VRmax-0.25.*(VRmax-VRmin).*rand(1,D));
                    
            e(i)=feval(fhd, pos(i,:)', varargin{:});% Evaluate fitness   
            fitcount=fitcount+1;

            if  e(i)<pbest_val(i) % update pbest value and position
                pbest(i,:)=pos(i,:);   
                pbest_val(i)=e(i);
                obj_func_slope(i)=0;
            else
                obj_func_slope(i)=obj_func_slope(i)+1;
            end
            if  pbest_val(i)<gbestval % update gbest value and postion
                gbest=pbest(i,:); 
                gbestval=pbest_val(i);
            end   
    
     end % CLPSO        
       
      for i=num_g1+1:num_g1+num_g2
                  
            if k <=1
                %pk = ones(1,5)*1/num_strategy;
                rk = 0:1/num_strategy:1;
                success_mem = zeros(1,num_strategy);
                failure_mem = zeros(1,num_strategy);
                %sk = zeros(1,num_strategy);
            elseif mod(fitcount,LP)==0
                total = (success_mem+failure_mem);
                total(total==0)=1;
                sk = (success_mem./total)+0.01;
                pk = sk./sum(sk);
                rk = [0 cumsum(pk)];
                success_mem = zeros(1,num_strategy);
                failure_mem = zeros(1,num_strategy);   
            end  
            probability = rand(1);
             
            if probability>= rk(1) && probability < rk(2)
                % PSO  
                strategy_k = 1;
                delta(i,:)=c2_1(k).*rand(1,D).*(pbest(i,:)-pos(i,:))+c2_2(k).*rand(1,D).*(gbestrep(i,:)-pos(i,:));
                vel(i,:)=w2(k).*vel(i,:)+delta(i,:); 
                
            elseif probability>= rk(2) && probability < rk(3)   
                % FDR-PSO
                strategy_k=2;
                dis=abs(repmat(pbest(i,:),num_g2,1)-pbest(1:num_g2,:));
                fiterr=repmat(pbest_val(i),num_g2,1)-pbest_val(1:num_g2)';
                fiterr=repmat(fiterr,1,D);
                fiterr=fiterr-(dis==zeros(num_g2,D)).*fiterr;
                dis=dis+(dis==zeros(num_g2,D));
                FDR=fiterr./dis;
                [~,Fid]=max(FDR);
                for dimcnt=1:D
                    Pnd(i,dimcnt)=pbest(Fid(dimcnt),dimcnt);
                end
                delta(i,:)=fii(1).*rand(1,D).*(pbest(i,:)-pos(i,:))+fii(2).*rand(1,D).*(gbestrep(i,:)-pos(i,:))+fii(3).*rand(1,D).*(Pnd(i,:)-pos(i,:));
                vel(i,:)=w2(k).*vel(i,:)+delta(i,:); 
   
                
            elseif probability>= rk(3) && probability < rk(4)   
                % HPSO
                strategy_k = 3;
                vel(i,:)=(c3_1(k).*rand(1,D).*(pbest(i,:)-pos(i,:)))+(c3_2(k).*rand(1,D).*(gbestrep(i,:)-pos(i,:)));
                for d=1:D
                    if vel(i,d)==0
                        if (rand(1)<0.5)
                            vel(i,d)=rand(1)*re_init_vel(k);
                        else
                            vel(i,d)=-rand(1)*re_init_vel(k);
                        end
                    end
                    vel(i,d)=sign(vel(i,d))*min(abs(vel(i,d)),VRmax);
                end
             
            elseif probability>= rk(4) && probability < rk(5)
                    %LIPS
                    strategy_k = 4;
                    EU_dist=dist(pos(i,:),pbest'); 
                    EU_dist(i)=max(EU_dist);
                    [~,min_index]=sort(EU_dist); 
                    fi=(4.1./nsize).*rand(nsize,D);
                    FIP=sum(fi.*pbest(min_index(1:nsize),:))./sum(fi);  
                    delta(i,:)=sum(fi).*(FIP-pos(i,:));
                    vel(i,:)=0.7298.*(vel(i,:)+delta(i,:));
%              
           elseif probability>= rk(5) && probability < rk(6)
                    %CLPSO
                    strategy_k = 5;                          
                    delta(i,:)=(c4_1(k).*rand(1,D).*(fri_best_pos(i,:)-pos(i,:)))+(c4_2(k).*rand(1,D).*(gbestrep(i,:)-pos(i,:)));
                    vel(i,:)=w2(k)*vel(i,:)+delta(i,:);
                    
                     if obj_func_slope(i)>5
                       fri_best(i,:)=i*ones(1,D);
                       friend1=(ceil(ps*rand(1,D)));
                       friend2=(ceil(ps*rand(1,D)));
                       friend=(pbest_val(friend1)<pbest_val(friend2)).*friend1+(pbest_val(friend1)>=pbest_val(friend2)).*friend2;
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
            end

            % for all
            vel(i,:)=((vel(i,:)<Vmin).*Vmin)+((vel(i,:)>Vmax).*Vmax)+(((vel(i,:)<Vmax)&(vel(i,:)>Vmin)).*vel(i,:));
            der(i,:)= alpha*der(i,:) + (gbestrep(i,:)-last_gbest(i,:));
            last_gbest(i,:) = gbestrep(i,:);
            pos(i,:)=pos(i,:)+vel(i,:)+kd*der(i,:); 
            
            pos(i, :) =   ((pos(i,:)<=VRmax)&(pos(i,:)>=VRmin)) .* pos(i,:) ...
                        + (pos(i,:)<VRmin) .* (VRmin+0.25.*(VRmax-VRmin).*rand(1,D)) ...
                        + (pos(i,:)>VRmax) .* (VRmax-0.25.*(VRmax-VRmin).*rand(1,D));
                    
            e(i)=feval(fhd, pos(i,:)', varargin{:});% Evaluate fitness   
            fitcount=fitcount+1;

            if  e(i)<pbest_val(i) % update pbest value and position
                pbest(i,:)=pos(i,:);   
                pbest_val(i)=e(i);
                success_mem(strategy_k) = success_mem(strategy_k) +1;
            else
                failure_mem(strategy_k) = failure_mem(strategy_k) + 1;
            end
            if strategy_k==5 && e(i)<pbest_val(i)
                obj_func_slope(i)=0;
            else
                obj_func_slope(i)=obj_func_slope(i)+1;
            end
            if  pbest_val(i)<gbestval % update gbest value and postion
                gbest=pbest(i,:); 
                gbestval=pbest_val(i);
            end   
    
      end

else
    
     for i=1:ps
                  
            if k <=1
                %pk = ones(1,5)*1/num_strategy;
                rk = 0:1/num_strategy:1;
                success_mem = zeros(1,num_strategy);
                failure_mem = zeros(1,num_strategy);
                %sk = zeros(1,num_strategy);
            elseif mod(fitcount,LP)==0
                total = (success_mem+failure_mem);
                total(total==0)=1;
                sk = (success_mem./total)+0.01;
                pk = sk./sum(sk);
                rk = [0 cumsum(pk)];
                success_mem = zeros(1,num_strategy);
                failure_mem = zeros(1,num_strategy);   
            end   

            probability = rand(1);
           
            if probability>= rk(1) && probability < rk(2)
                % PSO  
                strategy_k = 1;
                delta(i,:)=c2_1(k).*rand(1,D).*(pbest(i,:)-pos(i,:))+c2_2(k).*rand(1,D).*(gbestrep(i,:)-pos(i,:));
                vel(i,:)=w2(k).*vel(i,:)+delta(i,:);              

            elseif probability>= rk(2) && probability < rk(3)
                % FDR_PSO  
                strategy_k = 2;
                dis=abs(repmat(pbest(i,:),num_g2,1)-pbest(1:num_g2,:));
                fiterr=repmat(pbest_val(i),num_g2,1)-pbest_val(1:num_g2)';
                fiterr=repmat(fiterr,1,D);
                fiterr=fiterr-(dis==zeros(num_g2,D)).*fiterr;
                dis=dis+(dis==zeros(num_g2,D));
                FDR=fiterr./dis;
                [~,Fid]=max(FDR);
                for dimcnt=1:D
                    Pnd(i,dimcnt)=pbest(Fid(dimcnt),dimcnt);
                end
                delta(i,:)=fii(1).*rand(1,D).*(pbest(i,:)-pos(i,:))+fii(2).*rand(1,D).*(gbestrep(i,:)-pos(i,:))+fii(3).*rand(1,D).*(Pnd(i,:)-pos(i,:));
                vel(i,:)=w2(k).*vel(i,:)+delta(i,:); 
                
            elseif probability>= rk(3) && probability < rk(4)   
                % HPSO
                strategy_k = 3;
                vel(i,:)=(c3_1(k).*rand(1,D).*(pbest(i,:)-pos(i,:)))+(c3_2(k).*rand(1,D).*(gbestrep(i,:)-pos(i,:)));
                for d=1:D
                    if vel(i,d)==0
                        if (rand(1)<0.5)
                            vel(i,d)=rand(1)*re_init_vel(k);
                        else
                            vel(i,d)=-rand(1)*re_init_vel(k);
                        end
                    end
                    vel(i,d)=sign(vel(i,d))*min(abs(vel(i,d)),VRmax);
                end
             
            elseif probability>= rk(4) && probability < rk(5)
                    %LIPS
                    strategy_k = 4;
                    EU_dist=dist(pos(i,:),pbest'); 
                    EU_dist(i)=max(EU_dist);
                    [~,min_index]=sort(EU_dist); 
                    fi=(4.1./nsize).*rand(nsize,D);
                    FIP=sum(fi.*pbest(min_index(1:nsize),:))./sum(fi);  
                    delta(i,:)=sum(fi).*(FIP-pos(i,:));
                    vel(i,:)=0.7298.*(vel(i,:)+delta(i,:));
                    
           elseif probability>= rk(5) && probability < rk(6)
                    %CLPSO
                    strategy_k = 5; 
                    delta(i,:)=(c4_1(k).*rand(1,D).*(fri_best_pos(i,:)-pos(i,:)))+(c4_2(k).*rand(1,D).*(gbestrep(i,:)-pos(i,:)));
                    vel(i,:)=w2(k)*vel(i,:)+delta(i,:);
                                   
                    if obj_func_slope(i)>5
                       fri_best(i,:)=i*ones(1,D);
                       friend1=(ceil(ps*rand(1,D)));
                       friend2=(ceil(ps*rand(1,D)));
                       friend=(pbest_val(friend1)<pbest_val(friend2)).*friend1+(pbest_val(friend1)>=pbest_val(friend2)).*friend2;
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
            end
            
            % for all
            vel(i,:)=((vel(i,:)<Vmin).*Vmin)+((vel(i,:)>Vmax).*Vmax)+(((vel(i,:)<Vmax)&(vel(i,:)>Vmin)).*vel(i,:));
            der(i,:)= alpha*der(i,:) + (gbestrep(i,:)-last_gbest(i,:));
            last_gbest(i,:) = gbestrep(i,:);
            pos(i,:)=pos(i,:)+vel(i,:)+kd*der(i,:); 
            
            pos(i, :) =   ((pos(i,:)<=VRmax)&(pos(i,:)>=VRmin)) .* pos(i,:) ...
                        + (pos(i,:)<VRmin) .* (VRmin+0.25.*(VRmax-VRmin).*rand(1,D)) ...
                        + (pos(i,:)>VRmax) .* (VRmax-0.25.*(VRmax-VRmin).*rand(1,D));
                    
            e(i)=feval(fhd, pos(i,:)', varargin{:});% Evaluate fitness   
            fitcount=fitcount+1;

            if  e(i)<pbest_val(i) % update pbest value and position
                pbest(i,:)=pos(i,:);   
                pbest_val(i)=e(i);
                success_mem(strategy_k) = success_mem(strategy_k) +1;
            else
                failure_mem(strategy_k) = failure_mem(strategy_k) + 1;
            end
            if strategy_k==5 && e(i)<pbest_val(i)
                obj_func_slope(i)=0;
            else
                obj_func_slope(i)=obj_func_slope(i)+1;
            end
            if  pbest_val(i)<gbestval % update gbest value and postion
                gbest=pbest(i,:); 
                gbestval=pbest_val(i);
            end   

     end
end
hist_p(k+1,:) = gbest;
hist_v(k+1) = gbestval; 

end
FES=fitcount;
end