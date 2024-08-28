import matplotlib.pyplot as plt
import pandas as pd


data = {
    'P': [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
    'N': ['8192^2', '8192^2', '8192^2', '8192^2', '8192^2', '8192^2', '8192^2', '8192^2', '8192^2', '8192^2',
          '8192^2', '8192^2', '8192^2', '8192^2', '8192^2', '8192^2', '8192^2', '8192^2', '8192^2', '8192^2'],
    'Method': ['leader', 'noleader', 'leader', 'noleader', 'leader', 'noleader', 'leader', 'noleader', 'leader', 'noleader',
               'leader', 'noleader', 'leader', 'noleader', 'leader', 'noleader', 'leader', 'noleader', 'leader', 'noleader'],
    'Time': [19.006281, 19.172385, 19.398364, 24.748926, 21.978163, 25.665574, 16.841819, 19.880733, 24.487058, 26.740307,
             22.554482, 26.389119, 18.783649, 28.349476, 25.248739, 27.129253, 26.855486, 27.961464, 23.485604, 32.803622]
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
