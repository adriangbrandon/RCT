import matplotlib.pyplot as plt
import numpy as np

mpl_fig = plt.figure()
ax = mpl_fig.add_subplot(111)


ref_move = []
ref_rmq = []
traj_disap = []
traj_length = []
traj_offset = []
traj_values = []
traj_rmq = []
lists = [[6627002, 173820,43869631, 17535347, 124452576, 249043360, 7778286], 
		 [7929784, 304520,43869631, 16085014, 109286688, 218711584, 6830418], 
		 [9271699, 466043,43869631, 15237075, 100778336, 201694880, 6298646],
		 [10634254, 584826,43869631, 14794255, 96427776, 192993760, 6026736],
		 [11733475, 688073,43869631, 14646107, 95249856, 190637920, 5953116],
		 [13196673, 858373,43869631, 14235435, 91277440, 182693088, 5704840],
		 [14452537, 970264,43869631, 14077495, 89919776, 179977760, 5619986],
		 [15881522, 1091894,43869631, 13839233, 87621664, 175381536, 5476354],
		 [16990601, 1224930,43869631, 13771488, 87104064, 174346336, 5444004]
		 ]
for l in lists:
	ref_move.append(l[0]/8)
	ref_rmq.append(l[1]/8)
	traj_disap.append(l[2]/8)
	traj_length.append(l[3]/8)
	traj_offset.append(l[4]/8)
	traj_values.append(l[5]/8)
	traj_rmq.append(l[6]/8)

N = len(lists)
ind = np.arange(N)    # the x locations for the groups
width = 0.35       # the width of the bars: can also be len(x) sequence

p1 = ax.bar(ind, ref_move, width, color='red', label='ref_move')

bottom_array = ref_move;
p2 = ax.bar(ind, ref_rmq, width, color='blue', label='ref_rmq',
             bottom=bottom_array)

bottom_array = [sum(x) for x in zip(ref_rmq, bottom_array)]
p3 = ax.bar(ind, traj_disap, width, color='orange', label='traj_disap',
             bottom=bottom_array)

bottom_array = [sum(x) for x in zip(traj_disap, bottom_array)]
p4 = ax.bar(ind, traj_length, width, color='gray', label='traj_length',
             bottom=bottom_array)

bottom_array = [sum(x) for x in zip(traj_length, bottom_array)]
p5 = ax.bar(ind, traj_offset, width, color='brown', label='traj_offset',
             bottom=bottom_array)

bottom_array = [sum(x) for x in zip(traj_offset, bottom_array)]
p6 = ax.bar(ind, traj_values, width, color='green', label='traj_values',
             bottom=bottom_array)

bottom_array = [sum(x) for x in zip(traj_values, bottom_array)]
p7 = ax.bar(ind, traj_rmq, width, color='cyan', label='traj_rmq',
             bottom=bottom_array)

p8 = ax.bar([N+1],  [32153664], width, color='black')

ax.set_ylabel('Size (bytes)')
ax.set_xlabel('RCT reference')
ax.set_title('Block size = 4096')

ax.set_xticks(ind)
ax.set_yticks(np.arange(0, 90000000, 10000000))
ax.set_xticklabels(('1', '2', '3', '4', '5', '6', '7', '8', '9'))

#chartBox = ax.get_position()
#ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.6, chartBox.height])
ax.legend(loc='upper right', ncol=2)
plt.show();

#plotly_fig = tls.mpl_to_plotly( mpl_fig )

# For Legend
#plotly_fig["layout"]["showlegend"] = True
#plotly_fig["data"][0]["name"] = "Men"
#plotly_fig["data"][1]["name"] = "Women"
#py.iplot(plotly_fig, filename='stacked-bar-chart')