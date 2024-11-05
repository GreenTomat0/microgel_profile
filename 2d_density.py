import numpy as np
import math
import matplotlib.pyplot as plt
# from mpl_toolkits import mplot3d
from scipy.signal import savgol_filter



# (x0, y0) - center
x0 = 92#17#
y0 = 92#20#

x = []
y = []
gray = []


with open('central_particle_80th_slice.txt') as f: #    #small_square  #central_particle_80th_slice #sub1-5_xy_80th_slice_bg-401
	for row in f:
		x.append( float(row.split('\t')[0]) - x0 )
		y.append( float(row.split('\t')[1]) - y0 )
		gray.append( float(row.split('\t')[2]) )
f.close()



vectors_num = len(x)

r = [None]*vectors_num


for i in range(vectors_num):
	r[i] = math.sqrt(x[i]*x[i] + y[i]*y[i])



zipped_lists = zip(r, gray)
sorted_pairs = sorted(zipped_lists)

tuples = zip(*sorted_pairs)
r, gray = [ list(tuple) for tuple in  tuples]



r_new = []
gray_new = []

step = 1
eps = 1
if_not_end_list = True
i = 0

while if_not_end_list:
	r_new.append(r[i])
	gray_var = gray[i]

	count = 1
	j = i + 1
	while abs(r[i] - r[j]) < eps:
		count += 1
		gray_var += gray[j]
		j += 1

	gray_new.append(gray_var/count)
	
	i = j + 1

	if i >= len(r):
		if_not_end_list = False



indices_delete = [inx for inx, val in enumerate(r_new) if val > x0]
del r_new[indices_delete[0]:]
del gray_new[indices_delete[0]:]



#----- rearrange data -----#
# normilize to mingray value
min_gray_value = min(gray_new)
gray_new = [elem - min_gray_value for elem in gray_new]
# pixels to nm
multipl_coef = 470/184 # 470 nm (or 184 px) is the size of window
r_new = [elem*multipl_coef for elem in r_new]



# smooth data
gray_new_smooth = savgol_filter(gray_new, 51, 3)



f_2d = open("data/2d_data.dat", "w")
f_2d.write('r\tgray\tsmoothed\n')

for i in range(len(r_new)):
	f_2d.write(f'{r_new[i]}\t{gray_new[i]}\t{gray_new_smooth[i]}\n')
f_2d.close()



fig, axs = plt.subplots(figsize=(11, 7))


axs.plot( r_new, gray_new, c='gray', ls="none", marker='*', markerfacecolor='black', markersize=8, label='net data' )
axs.plot( r_new, gray_new_smooth, c='blue', lw=3, label='smoothed' )
# axs.plot( r_new, gray_new, c='black', ls="none", marker='o', markerfacecolor='gray', markersize=5, label='net data' )
# axs.plot( r_new, gray_new_smooth, c='red', lw=3, label='smoothed' )

plt.xlabel('r [nm]', fontsize=20)
plt.ylabel('relative gray intensity', fontsize=20) # mean gray value
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
axs.legend(fontsize=16, loc='best')

plt.title('2D', fontsize=22)

plt.savefig('images/mean-gray-value_80th_slice_2d.jpg', dpi=100)

plt.show()

# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.scatter(x, y, gray) 
# plt.show()
