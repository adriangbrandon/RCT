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
lists = [[6638105, 184793,43869631, 17147190, 119549216, 239236640, 7471826], 
		 [7969250, 316900,43869631, 15848076, 106404384, 212946976, 6650274], 
		 [9254033, 452304,43869631, 15216301, 100434624, 201007456, 6277164],
		 [10574424, 573651,43869631, 14793838, 96394880, 192927968, 6024680],
		 [11827781, 694901,43869631, 14519131, 93842688, 187823584, 5865168],
		 [13181003, 848142,43869631, 14226310, 91136864, 182411936, 5696054],
		 [14395997, 951126,43869631, 14084151, 89871808, 179881824, 5616988],
		 [15811098, 1101276,43869631, 13851701, 87744480, 175627168, 5484030],
		 [17010967, 1209633,43869631, 13745019, 86808416, 173755040, 5425526] 
		 ]
for l in lists:
	ref_move.append(l[0]/8)
	ref_rmq.append(l[1]/8)
	traj_disap.append(l[2]/8)
	traj_length.append(l[3]/8)
	traj_offset.append(l[4]/8) #*0.75
	traj_values.append(l[5]/8) #*0.5
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
ax.set_title('Block size = 1024')

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