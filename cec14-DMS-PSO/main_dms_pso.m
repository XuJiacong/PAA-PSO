clear;clc;
%mex cec14_func.cpp -DWINDOWS
func_num = 25;
D = 30;
Xmin = -100;
Xmax = 100;
pop_size= 20;
max_fes = 10000*D;
iter_max = ceil(max_fes/pop_size)-1100;
runs = 51;
fhd = str2func('cec14_func'); 
method = 1;
hist_v_all = zeros(iter_max+1, runs);
if method == 1
    for j=1:runs
            j,
            [gbest,gbestval,hist_v,hist_p,FES]= DMS_PSO_func(fhd,D,4,ceil(pop_size/4),iter_max,Xmin,Xmax,func_num);
            hist_v_all(:,j)=hist_v;
            xbest(j,:)=gbest;
            fbest(j)=gbestval;
            fbest(j)
    end
    f_mean=mean(fbest);
    f_std = std(fbest);
    hist_v = median(hist_v_all, 2);
elseif method == 2
    for j=1:runs
            j,
            [gbest,gbestval,hist_v,hist_p,FES]= PBS_DMS_PSO_func(fhd,D,4,ceil(pop_size/4),iter_max,Xmin,Xmax,func_num);
            hist_v_all(:,j)=hist_v;
            xbest(j,:)=gbest;
            fbest(j)=gbestval;
            fbest(j)
    end
    f_mean = mean(fbest);
    f_std = std(fbest);
    hist_v = median(hist_v_all, 2);
elseif method == 3
    for j=1:runs
            j,
            [gbest,gbestval,hist_v,hist_p,FES]= PID_DMS_PSO_funcM(fhd,D,4,ceil(pop_size/4),iter_max,Xmin,Xmax,func_num);
            hist_v_all(:,j)=hist_v;
            xbest(j,:)=gbest;
            fbest(j)=gbestval;
            fbest(j)
    end
    f_mean=mean(fbest);
    f_std = std(fbest);
    hist_v = median(hist_v_all, 2);
elseif method == 4
    for j=1:runs
            j,
            [gbest,gbestval,hist_v,hist_p,FES]= PID_DMS_PSO_func(fhd,D,4,ceil(pop_size/4),iter_max,Xmin,Xmax,func_num);
            hist_v_all(:,j)=hist_v;
            xbest(j,:)=gbest;
            fbest(j)=gbestval;
            fbest(j)
    end
    f_mean=mean(fbest);
    f_std = std(fbest);
    hist_v = median(hist_v_all, 2);

end
