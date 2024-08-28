import matplotlib.pyplot as plt

data = [
    [(512, '5 Point Stencil', 0.074896), (512, '5 Point Stencil', 0.064020), (512, '5 Point Stencil', 0.072064)],
    [(512, '9 Point Stencil', 0.082343), (512, '9 Point Stencil', 0.079425), (512, '9 Point Stencil', 0.083919)],
]

times_5_point_stencil = [time for _, stencil, time in data[0] if stencil == '5 Point Stencil']
times_9_point_stencil = [time for _, stencil, time in data[1] if stencil == '9 Point Stencil']

plt.boxplot([times_5_point_stencil, times_9_point_stencil], labels=['512, 5 Point Stencil', '512, 9 Point Stencil'])
plt.title('Execution Time for D=512')
plt.xlabel('Stencil Configuration')
plt.ylabel('Time (seconds)')

for i, times in enumerate([times_5_point_stencil, times_9_point_stencil]):
    for j, time in enumerate(times):
        plt.text(i + 1.1, time, f'{time:.6f} s', ha='left', va='center', fontsize=6)


plt.show()
