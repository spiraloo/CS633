import matplotlib.pyplot as plt
import pandas as pd


data = {
    'P': [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
    'N': ['4096^2', '4096^2', '4096^2', '4096^2', '4096^2', '4096^2', '4096^2', '4096^2', '4096^2', '4096^2',
          '4096^2', '4096^2', '4096^2', '4096^2', '4096^2', '4096^2', '4096^2', '4096^2', '4096^2', '4096^2'],
    'Method': ['leader', 'noleader', 'leader', 'noleader', 'leader', 'noleader', 'leader', 'noleader', 'leader', 'noleader',
               'leader', 'noleader', 'leader', 'noleader', 'leader', 'noleader', 'leader', 'noleader', 'leader', 'noleader'],
    'Time': [4.882410, 4.924945, 5.011626, 5.031206, 4.690415, 4.759920, 4.873025, 5.339797, 5.455722, 5.490392,
             4.980741, 5.065989, 4.873025, 5.339797, 6.030327, 6.147957, 5.754658, 5.804026, 5.459630, 5.475138]
}


df = pd.DataFrame(data)

# Plot
ax = df.boxplot(column='Time', by=['P', 'N', 'Method'], grid=False, figsize=(12, 8))
ax.set_title('Time Taken for Each Data Size per Method per Process Count')
ax.set_xlabel('(P, N, Method)')
ax.set_ylabel('Time (seconds)')
plt.xticks(rotation=30)


plt.subplots_adjust(bottom=0.2)

plt.show()
