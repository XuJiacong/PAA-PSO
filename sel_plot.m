clear;clc;
std1 = "GPSO";
root1 = "cec14-SPSO/";
std2 = "HCLPSO";
root2 = "cec14-HCLPSO/";
me = 15000;
sel_set = [2,4,6,7,8,9,10,11,13,14,15,18,23,24,29,30];
gap = 500;
for fc = sel_set
    
    %%%% GPSO %%%%
    filename = root1+"func"+fc+"_org_"+std1+".mat";
    load(filename)
    plot(log(hist_v),'-o','MarkerIndices',1:gap:length(hist_v),'Color', 'r', 'LineWidth', 1); hold on;
    
    
    filename = root1+"func"+fc+"_PBS_"+std1+".mat";
    load(filename)
    plot(log(hist_v),'-d','MarkerIndices',1:gap:length(hist_v),'Color', 'r', 'LineWidth', 1); hold on;

    filename = root1+"func"+fc+"_PIDM_"+std1+".mat";
    load(filename)
    plot(log(hist_v),'-+','MarkerIndices',1:gap:length(hist_v),'Color', 'r', 'LineWidth', 1); hold on;

    filename = root1+"func"+fc+"_PID_"+std1+".mat";
    load(filename)
    plot(log(hist_v),'-^','MarkerIndices',1:gap:length(hist_v),'Color', 'r', 'LineWidth', 1); hold on;
    
    %%%% HCLPSO %%%%
    filename = root2+"func"+fc+"_org_"+std2+".mat";
    load(filename)
    plot(log(hist_v),'-.o','MarkerIndices',1:gap:length(hist_v),'Color', 'b', 'LineWidth', 1); hold on;
    
    
    filename = root2+"func"+fc+"_PBS_"+std2+".mat";
    load(filename)
    plot(log(hist_v),'-.d','MarkerIndices',1:gap:length(hist_v),'Color', 'b', 'LineWidth', 1); hold on;

    filename = root2+"func"+fc+"_PIDM_"+std2+".mat";
    load(filename)
    plot(log(hist_v),'-.+','MarkerIndices',1:gap:length(hist_v),'Color', 'b', 'LineWidth', 1); hold on;

    filename = root2+"func"+fc+"_PID_"+std2+".mat";
    load(filename)
    plot(log(hist_v),'-.^','MarkerIndices',1:gap:length(hist_v),'Color', 'b', 'LineWidth', 1); hold on;

    %title("Function "+fc, 'Interpreter','latex')
    legend(std1,"PBS-"+std1,"PBSv2-"+std1,"PAA-"+std1,std2,"PBS-"+std2,"PBSv2-"+std2,"PAA-"+std2)
    xlabel('Iteration', 'Interpreter','latex')
    ylabel('log(Median Fitness)','Interpreter','latex')
    
    hold off;
    
    savedir = "figures/cec14-Mixed/"+"func"+fc+"_all"+".jpg";
    saveas(gcf, savedir)
    
end