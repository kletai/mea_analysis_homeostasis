function chosen_index = sel_index_popup()
    prompt = {'Enter electrode index:'};
    dlg_title = 'Choose Electrode';
    num_lines = 1;
    chosen_index = inputdlg(prompt, dlg_title, num_lines);
    if isempty(chosen_index)
        chosen_index = 0;
    else
        chosen_index = str2num(chosen_index{1});
    end
    disp(['chose = ' num2str(chosen_index)]);