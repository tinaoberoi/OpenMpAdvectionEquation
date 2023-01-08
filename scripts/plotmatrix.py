import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
import os

file_name1 = 'timestamp_8000.txt'
#file_name2 = 'timestamp2.txt'
path1 = os.path.abspath(file_name1)
#path2 = os.path.abspath(file_name2)
data_time1 = np.genfromtxt(path1, delimiter=',')
#data_time2 = np.genfromtxt(path2, delimiter=',')
plt.imshow(data_time1, interpolation='none')
plt.savefig('timestamp_8000.png')
#plt.imshow(data_time2, interpolation='none')
#plt.savefig('test2.png')
#plt.show()
print("something")