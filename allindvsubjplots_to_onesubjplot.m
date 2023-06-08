%function [] = allindvsubjplots_to_onesubjplot(subjidx, fitted_on_all_data)
        subjidx=7; fitted_on_all_data=false; fontsize = 9;
        
   
    if(~fitted_on_all_data)
        figure(1)
        nexttile(subjidx);
        ax1=gca;
        figure(2)
        nexttile(subjidx);
        ax2=gca;

        figure;
        T=tiledlayout(1,2,'Padding', 'compact', 'TileSpacing', 'compact');
        t1 = nexttile(1);
        hold on
        fig1 = get(ax1,'children');
        copyobj(fig1, t1);
        xlabel("Stimulus location (\circ)", 'FontSize', fontsize)
        ylabel("Mean loc estimate (\circ)", 'FontSize', fontsize)
        ylim([-20,20])
        lg = legend("Visual (high reliability)","Visual (med reliability)", "Visual (low reliability)", "Auditory");
        set(lg,'Box','off')
        
        t2 = nexttile(2);
        fig2 = get(ax2,'children');
        copyobj(fig2, t2);
        xlabel("Stimulus location (\circ)", 'FontSize', fontsize)
        ylabel("SD of loc estimate (\circ)", 'FontSize', fontsize)
        ylim([-20,20])
        lg = legend("Visual (high reliability)","Visual (med reliability)", "Visual (low reliability)", "Auditory");
        set(lg,'Box','off')
        
        
        
    else
    end
%end