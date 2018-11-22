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
lists = [[6635213, 188158,43869631, 17091627, 118602304, 237342816, 7412644], 
		 [7942544, 322989,43869631, 15873429, 106509600, 213157408, 6656850], 
		 [9233798, 450427,43869631, 15259179, 100639904, 201418016, 6289994],
		 [10539551, 581888,43869631, 14849592, 96797184, 193732576, 6049824],
		 [11829496, 720790,43869631, 14527274, 93827104, 187792416, 5864194]
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
ax.set_title('Block size = 256')

ax.set_xticks(ind + width/2.)
ax.set_yticks(np.arange(0, 90000000, 10000000))
ax.set_xticklabels(('1', '2', '3', '4', '5'))

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