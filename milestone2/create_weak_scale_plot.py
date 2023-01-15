import matplotlib.pyplot as plt
import numpy as np

# for  200 N guided lax order
# first point is serial
xpoints = np.array([1, 2, 3, 4, 5, 6, 7])
ypoints = np.array([3, 8, 1, 10])
# serial, (2), (4), (6), (8), (10), (12), (16), (20), (32), (45), (64)
# first_order
time = np.array([
    0.281484,
    0.287624,
    0.286564,
    0.293673,
    0.348915,
    0.980967,
    1.701600
])
# lax
# time = np.array([55.132473/55.132473,

# 55.132473/27.812475,

# 55.132473/13.910441,

# 55.132473/9.366285,

# 55.132473/7.131303,

# 55.132473/5.673626,

# 55.132473/6.499164,

# 55.132473/5.477798,

# 55.132473/6.562273,

# 55.132473/6.407046,

# 55.132473/6.243339,

# 55.132473/6.686601

# ])
# second_order
#time = np.array([143.761289/143.761289, 71.385954/143.761289, 35.537750/143.761289, 23.799928/143.761289, 17.965516/143.761289, 11.958249/143.761289, 9.291335/143.761289, 10.708042/143.761289, 10.610915/143.761289, 11.369717/143.761289, 10.907213/143.761289 ])
# guassian
#time = np.array([0.184164, 0.096234/0.184164, 0.050016/0.184164, 0.037531/0.184164, 0.029898/0.184164, 0.025289/0.184164, 0.028791/0.184164, 0.018069/0.184164, 0.019931/0.184164, 0.023091/0.184164, 0.025380/0.184164])

labels = ['N = 200 n = 1', 'N = 283 n = 2', 'N = 400 n = 4', 'N = 565 n = 8', 'N = 800 n = 16', 'N = 1131 n = 32', 'N = 1600 n = 64']
plt.xlabel("Number Of Threads")
plt.ylabel('Time Taken')
plt.title("Weak Scale Static Scheduling Performance")
plt.text(-5, 60, 'N = Grid Size', fontsize = 12)
plt.text(5, 60, 'n = Number of cores', fontsize = 12)
plt.xticks(xpoints, labels, rotation ='vertical')
plt.plot(xpoints, time)
plt.savefig('weak_scale.png')
plt.show()