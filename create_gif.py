import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os
nframes = 10
plt.subplots_adjust(top=1, bottom=0, left=0, right=1)

def animate(i):
    file_name1 = 'plots/' + 'image' + str(i) + '.png'
    #file_name2 = 'timestamp2.txt'
    path1 = os.path.abspath(file_name1)
    im = plt.imread(path1)
    plt.imshow(im)

anim = FuncAnimation(plt.gcf(), animate, frames=nframes, interval=500)
anim.save('output.gif', writer='imagemagick')