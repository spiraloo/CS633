import matplotlib.pyplot as plt


data = [
    [(2048, '5 Point Stencil', 0.859955), (2048, '5 Point Stencil', 0.756573), (2048, '5 Point Stencil', 1.092581)],
    [(2048, '9 Point Stencil', 1.254523), (2048, '9 Point Stencil',0.882341), (2048, '9 Point Stencil',1.237472)],
]


times_5_point_stencil = [time for _, stencil, time in data[0] if stencil == '5 Point Stencil']
times_9_point_stencil = [time for _, stencil, time in data[1] if stencil == '9 Point Stencil']


plt.boxplot([times_5_point_stencil, times_9_point_stencil], labels=['5 Point Stencil', '9 Point Stencil'])
plt.title('Execution Time for D=2048')
plt.xlabel('Stencil Configuration')
plt.ylabel('Time (seconds)')


for i, times in enumerate([times_5_point_stencil, times_9_point_stencil]):
    for j, time in enumerate(times):
        plt.text(i + 1.1, time, f'{time:.6f} s', ha='left', va='center', fontsize=6)

plt.show()