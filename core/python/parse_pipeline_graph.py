import json

def uniquify_list(seq):  
    checked = [] 
    for e in seq: 
        if e not in checked: 
            checked.append(e) 
    return checked

def frame_type_str_to_color(s):
    cmap = {'Timepoint': 'm',
            'Housekeeping' : 'y',
            'Scan' : 'r',
            'Map' : 'g',
            'InstrumentStatus' : 'k',
            'Wiring' : 'c',
            'Calibration' : 'b',
            'EndProcessing' : 'k',
            'None' : 'k',
        }
    assert(s in cmap)
    return cmap[s]

def plot_frame_processing_info(g3_pipeline):
    import matplotlib.pyplot as plt

    plt.clf()
    parsed_output = json.loads( g3_pipeline.GetGraphInfo()  )

    proc_lst_trans = zip(*parsed_output['Processing List'])
    ptime = proc_lst_trans[0]
    pmod = proc_lst_trans[1]
    pframe = proc_lst_trans[2]
    ptype = proc_lst_trans[3]

    n_frames = max(pframe)+1
    n_times = max(ptime)+1
    n_mods = max(pmod)+1


    #create the legend in a realllly hacky way
    unique_ptypes = uniquify_list(ptype)
    ax = plt.gca()
    colors = map(frame_type_str_to_color, unique_ptypes)
    labels = unique_ptypes
    [plt.plot([], [], color=c,label=l, linewidth=3 )[0] for c, l in zip(colors,labels)]
    plt.legend(labels)

    #these store the information for the state of the frame
    frame_mod_pos = {}
    frame_prev_time = {}

    #plots the data flow
    for i in range(len(ptime)):
        if pframe[i] not in frame_mod_pos:
            frame_mod_pos[pframe[i]] = pmod[i] - 1
        prev_mod = frame_mod_pos[pframe[i]]
        if pframe[i] in frame_prev_time:
            prev_time = frame_prev_time[pframe[i]]
        else:
            prev_time = ptime[i]
        cur_mod = pmod[i]
        cur_time = ptime[i]
        frame_id = pframe[i]
        x_pos = lambda mod: mod - float(frame_id + 1)/(3*float(n_frames + 1)) + .66
        plt.plot([x_pos(prev_mod), x_pos(prev_mod)], 
                 [prev_time, cur_time], 
                 color = frame_type_str_to_color(ptype[i]),
                 linestyle = '--'
        )
        plt.plot([x_pos(prev_mod), x_pos(cur_mod)], 
                 [cur_time, cur_time+1],
                 color = frame_type_str_to_color(ptype[i]),
                 linewidth = 2.0,
                 marker = '.'
             )
        frame_mod_pos[pframe[i]] = cur_mod
        frame_prev_time[pframe[i]] = cur_time + 1

    #handles setting up labels
    mod_x_ticks = []
    mod_x_labels = []

    module_stripper = lambda mods_name: mods_name[mods_name.rfind('.')+1:]
    for p in parsed_output["Module List"]:
        plt.axvline(p[0], linewidth = 2, linestyle = '-.', color = 'k')
        mod_x_ticks.append(p[0]+.5)
        mod_x_labels.append(module_stripper(p[1]))

    plt.xticks(mod_x_ticks, mod_x_labels, rotation=-7)    
    #does some silly formatting
    plt.gca().get_yaxis().set_ticks([])
    plt.ylim(n_times+.5, -0.5)
    plt.xlim(0, n_mods)
    plt.ylabel("Time\n$\\Longleftarrow$")
    plt.xlabel("Module $\\Longrightarrow$")
    plt.show()

