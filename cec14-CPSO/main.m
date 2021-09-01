clear;clc;  
%mex cec14_func.cpp -DWINDOWS
D = 30;
Xmin = -100;
Xmax = 100; 
pop_size= 20;
max_fes = 10000*D;
iter_max = ceil((max_fes-pop_size)/pop_size/7);
runs = 51;
fhd = str2func('cec14_func'); 
method = 4;
hist_v_all = zeros(iter_max, runs);
record = zeros(60, 1);
for i = 1:30
    func_num = i;
    if method == 1
        for j=1:runs
                j,
                [gbest,gbestval,hist_v,hist_p,FES]= CPSO_func(fhd,D,pop_size,iter_max,Xmin,Xmax,func_num);
                hist_v_all(:,j)=hist_v;
                xbest(j,:)=gbest;
                fbest(j)=gbestval;
                fbest(j)
        end
        f_mean=mean(fbest);
        f_std = std(fbest);
        hist_v = median(hist_v_all, 2);
        name = "func"+func_num+"_org_CPSO.mat";
        save(name, "hist_v")
        record(2*i-1) = f_mean;
        record(2*i) = f_std;
    elseif method == 2
        for j=1:runs
                j,
                [gbest,gbestval,hist_v,hist_p,FES]= PBS_CPSO_func(fhd,D,pop_size,iter_max,Xmin,Xmax,func_num);
                hist_v_all(:,j)=hist_v;
                xbest(j,:)=gbest;
                fbest(j)=gbestval;
                fbest(j)
        end
        f_mean = mean(fbest);
        f_std = std(fbest);
        hist_v = median(hist_v_all, 2);
        name = "func"+func_num+"_PBS_CPSO.mat";
        save(name, "hist_v")
        record(2*i-1) = f_mean;
        record(2*i) = f_std;
    elseif method == 3
        for j=1:runs
                j,
                [gbest,gbestval,hist_v,hist_p,FES]= PID_CPSO_funcM(fhd,D,pop_size,iter_max,Xmin,Xmax,func_num);
                hist_v_all(:,j)=hist_v;
                xbest(j,:)=gbest;
                fbest(j)=gbestval;
                fbest(j)
        end
        f_mean=mean(fbest);
        f_std = std(fbest);
        hist_v = median(hist_v_all, 2);
        name = "func"+func_num+"_PIDM_CPSO.mat";
        save(name, "hist_v")
        record(2*i-1) = f_mean;
        record(2*i) = f_std;
    else  
        for j=1:runs
                j,
                [gbest,gbestval,hist_v,hist_p,FES]= PID_CPSO_func(fhd,D,pop_size,iter_max,Xmin,Xmax,func_num);
                hist_v_all(:,j)=hist_v;
                xbest(j,:)=gbest;
                fbest(j)=gbestval;
                fbest(j)
        end
        f_mean=mean(fbest);
        f_std = std(fbest);
        hist_v = median(hist_v_all, 2);
        name = "func"+func_num+"_PID_CPSO.mat";
        save(name, "hist_v")
        record(2*i-1) = f_mean;
        record(2*i) = f_std;
    end
end
xlswrite('record.xlsx', record);