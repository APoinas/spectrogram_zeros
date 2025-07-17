import cv2
import os
os.chdir("/home/apoinas/Bureau/spectrogram_zeros")
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from Signal import Signal

def Remove_axis(ax):
    ax.xaxis.set_tick_params(labelbottom=False)
    ax.yaxis.set_tick_params(labelleft=False)
    ax.set_xticks([])
    ax.set_yticks([])

os.chdir("/home/apoinas/Bureau") # Folder in which to save the video
N_frame = 500 # Number of frames in the video
fps = 20 # FPS of the video
SNR = [200, 40] # SNR of both signals
b = 0.4 # Slope parameter of both signals
a_min = 0.05
a_max = 3 # Range for the distance between signals

a_list = np.linspace(-a_max, -a_min, N_frame)
video = cv2.VideoWriter('Anim.avi', cv2.VideoWriter_fourcc(*'DIVX'), fps, (720, 720))

for i in tqdm(range(N_frame)):
    np.random.seed(42)
    sig = Signal({"SNR":SNR, "a":[a_list[i], 0], "b":b}, T_length=15)
    fig, ax = sig.plot(show_max=True, plot_type=["spectro"])
    Remove_axis(ax[0])
    fig.canvas.draw()
    img_plot = np.array(fig.canvas.renderer.buffer_rgba())
    video.write(cv2.cvtColor(img_plot, cv2.COLOR_RGBA2BGR))
    plt.close()
video.release()
cv2.destroyAllWindows()