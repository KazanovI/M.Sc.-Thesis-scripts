function [shrd] = shared_sound_lick(ind_sound,ind_task)
    %% compare num of sig cells in start and task - and enh vs sup, per mouse
    inds = [ind_sound; ind_task];
    shrd = cell(1,length(inds));
    for s_i = 1:length(inds)
        shrd{s_i} = inds{1,s_i}' == inds{2,s_i} & inds{1,s_i}' == 1;
    end
end