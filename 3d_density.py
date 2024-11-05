import numpy as np
import math
import matplotlib.pyplot as plt
# from mpl_toolkits import mplot3d
from scipy.signal import savgol_filter


Nframes = 134 #92 # number of slices = z

# (x0, y0) - center
x0 = 92 #629.97/2 #
y0 = 92 #x0 #
# z0 = round(Nframes/2)


x = [] 
y = []
gray = []
#z = [i for i in range(-z0+1, z0)] #[i for i in range(-z0+1, z0)]

# print(z)


with open('sub1-5_xy_coord_134frames2.txt') as f:
	prev_row = -1 # index of previous row

	for row in f:
		cur_row = int(row.split('\t')[2]) # index of current row

		if cur_row != prev_row:
			x.append([])
			y.append([])
			gray.append([])
			
			prev_row = cur_row
		
		x[cur_row].append( float(row.split('\t')[0]) - x0 )
		y[cur_row].append( float(row.split('\t')[1]) - y0 )
		gray[cur_row].append( float(row.split('\t')[3]) )

f.close()


# print('rows:    ' + str(len(x)), '\n\n')

# for i in range(len(x)):
# 	print('columns:    ' + str(len(x[i])))

r = []

for i in range(Nframes):
	r.append([])
	for j in range(len(x[i])):
		
		r[i].append( math.sqrt(x[i][j]*x[i][j] + y[i][j]*y[i][j]) )


# print('rows:    ' + str(len(r)), '\n\n')

# for i in range(len(r)):
# 	print('columns:    ' + str(len(r[i])))


for i in range(Nframes):

	zipped_lists = zip(r[i], gray[i])
	sorted_pairs = sorted(zipped_lists)

	tuples = zip(*sorted_pairs)
	r[i], gray[i] = [ list(tuple) for tuple in  tuples]



# for i in [1, 23, 58, 94, 116]:
# 	for j in range(10):
# 		print(gray[i][j])
# 	print('\n\n')


r_new = []
gray_new = []

eps = 1

for k in range(Nframes):

	r_new.append([])
	gray_new.append([])

	# step = 1
	if_not_end_list = True
	i = 0

	while if_not_end_list:
		r_new[k].append(r[k][i])
		gray_var = gray[k][i]

		count = 1
		j = i + 1

		while abs(r[k][i] - r[k][j]) < eps:
			count += 1
			gray_var += gray[k][j]
			j += 1


		gray_new[k].append(gray_var/count)
		
		i = j + 1

		if i >= len(r[k]):
			if_not_end_list = False


	indices_delete = [inx for inx, val in enumerate(r_new[k]) if val > x0]
	# if len(indices_delete) > 0:
	del r_new[k][indices_delete[0]:]
	del gray_new[k][indices_delete[0]:]







max_lengths = []
max_length = 0

for i in range(Nframes):
	max_lengths.append(len(r_new[i]))
	if max_lengths[i] > max_length:
		max_length = max_lengths[i]



inx_delete = []

for i in range(Nframes):
	if len(r_new[i]) < max_length:
		inx_delete.append(i)



for i in sorted(inx_delete, reverse = True):
    del r_new[i]
    del gray_new[i]






# for i in range(len(r_new)):
# 	print(f'columns r_new[{i}]:   {len(r_new[i])}')

# print('\n\n')

# for i in range(len(gray_new)):
# 	print(f'columns gray_new[{i}]:   {len(gray_new[i])}')








r_3d = []
gray_3d = []


for k in range(len(r_new[0])):

	if_not_end_list = True
	i = 0
	while if_not_end_list:
		r_3d.append(r_new[i][k])
		gray_var = gray_new[i][k]

		count = 1
		j = i + 1
		
		while j < len(r_new) and abs(r_new[i][k] - r_new[j][k]) < eps:
			count += 1
			gray_var += gray_new[j][k]
			j += 1

		gray_3d.append(gray_var/count)
		
		i = j + 1

		if i >= len(r_new):
			if_not_end_list = False


print(len(r_3d))
print(len(gray_3d))



#----- rearrange data -----#
# normilize to mingray value
min_gray_value = min(gray_3d)
gray_3d = [elem - min_gray_value for elem in gray_3d]
# pixels to nm
multipl_coef = 470/184 # 470 nm (or 184 px) is the size of window
r_3d = [elem*multipl_coef for elem in r_3d]



gray_3d_smooth = savgol_filter(gray_3d, 51, 3)



r_2d = []
gray_2d = []
gray_2d_smooth = []
#----- read 2D data -----#
with open('data/2d_data.dat') as f_2d:
	next(f_2d)
	for row in f_2d:
		r_2d.append( float(row.split('\t')[0]) )
		gray_2d.append( float(row.split('\t')[1]) )
		gray_2d_smooth.append( float(row.split('\t')[2]) )
f_2d.close()



fig, axs = plt.subplots(figsize=(11, 7))

axs.plot( r_2d, gray_2d, c='gray', ls="none", marker='*', markerfacecolor='black', markersize=8, label='2D net data' )
axs.plot( r_2d, gray_2d_smooth, c='blue', lw=3, label='2D smoothed' )

axs.plot( r_3d, gray_3d, c='black', ls="none", marker='o', markerfacecolor='gray', markersize=5, label='3D net data' )
axs.plot( r_3d, gray_3d_smooth, c='red', lw=3, label='3D smoothed' )


plt.xlabel('r [nm]', fontsize=20)
plt.ylabel('relative gray intensity', fontsize=20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
axs.legend(fontsize=16, loc='best')

# plt.title('3D', fontsize=22)
plt.title('2D & 3D', fontsize=22)

# plt.savefig('images/mean-gray-value_134slides_3d.jpg', dpi=100)
# plt.savefig('images/center_particle_2d_3d.jpg', dpi=100)

plt.show()