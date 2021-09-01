clear;clc;
std = "GPSO";
root = "cec14-SPSO/";
me = 15000;
for fc = 1:30
    
    filename = root+"func"+fc+"_org_"+std+".mat";
    load(filename)
    plot(log(hist_v),'LineWidth', 2); hold on;
    h = plot(1:400:me,log(hist_v(1:400:me)),'O'); hold on;
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    filename = root+"func"+fc+"_PBS_"+std+".mat";
    load(filename)
    plot(log(hist_v),'LineWidth', 2); hold on;
    h = plot(1:400:me,log(hist_v(1:400:me)),'d'); hold on;
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';

    filename = root+"func"+fc+"_PIDM_"+std+".mat";
    load(filename)
    plot(log(hist_v),'LineWidth', 2); hold on;
    h = plot(1:400:me,log(hist_v(1:400:me)),'+'); hold on;
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';

    filename = root+"func"+fc+"_PID_"+std+".mat";
    load(filename)
    plot(log(hist_v),'LineWidth', 2); hold on;
    h = plot(1:400:me,log(hist_v(1:400:me)),'^'); hold on;
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';

    title("Function "+fc)
    legend(std,"PBS-"+std+"-V1","PBS-"+std+"-V2","PE-"+std)
    xlabel('Iteration')
    ylabel('log(Loss Value)')
    
    hold off;
    
    savedir = "figures/"+root+"func"+fc+"_all_"+std+".jpg";
    saveas(gcf, savedir)
    
end