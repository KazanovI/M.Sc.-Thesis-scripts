function plot_lick_snd_shr(ind_sound_lick_shared)

    lick_snd_shr = cell2mat(cellfun(@(x) sum(x),ind_sound_lick_shared,'UniformOutput',false)');
    perc_shr = (sum(lick_snd_shr)/sum(lick_snd_shr,'all'))*100;
    figure();
    bar(perc_shr,'FaceColor','k');
    xticklabels(["Enhanced","Suppressed","None"]);
    ylim([0 100]);
    ylabel('% of cells');
    box off
end