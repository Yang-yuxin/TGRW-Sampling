import matplotlib.pyplot as plt
import numpy as np

# Sample data
categories = ['1e4 walks', '1e5 walks', '1e6 walks']
values1 = np.log(1000*[0.0164, 0.163, 1.512])
values2 = np.log(1000*[0.0166, 0.161, 1.490])

# Create an array with the position of each bar along the x-axis
x_pos = np.arange(len(categories))

# Size of the bar
bar_width = 0.2

# Plotting
plt.bar(x_pos - bar_width/2, values1, width=bar_width, label='Runtime ITS')
plt.bar(x_pos + bar_width/2, values2, width=bar_width, label='Precomputed ITS')

# Adding labels and title
# plt.xlabel('')
plt.ylabel('log(Runtime) (ms)')
plt.title('Runtime of ITS and Precomputed ITS')
plt.xticks(x_pos, categories)

# Adding a legend
plt.legend()

# Display the plot
plt.savefig('pre_rt.png')
