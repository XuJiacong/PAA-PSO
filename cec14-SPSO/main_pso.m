clear;clc;
%mex cec14_func.cpp -DWINDOWS
func_num = 2;
D = 30;
Xmin = -100;
Xmax = 100;
pop_size= 20;
max_fes = 10000*D;
iter_max = int32(max_fes/pop_size);
runs = 51;
fhd = str2func('cec14_func'); 
method = 3;
hist_v_all = zeros(iter_max, 1);
if method == 1
    for j=1:runs
            j,
            [gbest,gbestval,hist_v,hist_p,FES]= PSO_func(fhd,D,pop_size,iter_max,Xmin,Xmax,func_num);
            hist_v_all = hist_v_all + hist_v./runs;
            xbest(j,:)=gbest;
            fbest(j)=gbestval;
            fbest(j)
    end
    f_mean=mean(fbest);
    f_std = std(fbest);
elseif method == 2
    for j=1:runs
            j,
            [gbest,gbestval,hist_v,hist_p,FES]= PBS_PSO_func(fhd,D,pop_size,iter_max,Xmin,Xmax,func_num);
            hist_v_all = hist_v_all + hist_v./runs;
            xbest(j,:)=gbest;
            fbest(j)=gbestval;
            fbest(j)
    end
    f_mean=mean(fbest);
    f_std = std(fbest);
elseif method == 3
    for j=1:runs
            j,
            [gbest,gbestval,hist_v,hist_p,FES]= PID_PSO_funcM(fhd,D,pop_size,iter_max,Xmin,Xmax,func_num);
            hist_v_all = hist_v_all + hist_v./runs;
            xbest(j,:)=gbest;
            fbest(j)=gbestval;
            fbest(j)
    end
    f_mean=mean(fbest);
    f_std = std(fbest);
else
    for j=1:runs
            j,
            [gbest,gbestval,hist_v,hist_p,FES]= PID_PSO_func(fhd,D,pop_size,iter_max,Xmin,Xmax,func_num);
            hist_v_all = hist_v_all + hist_v./runs;
            xbest(j,:)=gbest;
            fbest(j)=gbestval;
            fbest(j)
    end
    f_mean=mean(fbest);
end
