Python code accompanying the paper `On two fundamental properties of the zeros of
spectrograms of noisy signals <??????????????>`_ with 
`Rémi Bardenet <https://rbardenet.github.io/>`_.

The file "Signal.py" contains the "Signal" class with the function needed to simulate the signals and associated spectrograms. An exemple of utilisation of these functions is given in the file "Figures in the paper.ipynb", recreating the figures apppearing in the associated paper. The file "Generating_animations_video.py" contains some additional code to generate small videos showing how the zeros and local maxima of a noisy pair of paralle chirps behave when both chirps gets closer and closer from each other. An example of such a video is shown in the file "Anim.avi".

The Signal class
------------------------------------

Syntax
~~~~~~

The full syntax for defining an object of the Signal class is

Signal(par, win_length=512, T0=None, T_length=16, fs=1)

The main argument is "par" corresponding to a dictionnary of the parameters for the intended signal:

- For Hermite functions, "par" should be {"SNR":x, "k":y} where x is a positive float corresponding to the SNR of the signal and y is a positive int corresponding to the parameter of the Hermite function.
- For a single linear chirp, "par" should be {"SNR":x, "a":y, "b":z} where x is a positive float corresponding to the SNR of the signal and y and z are floats corresponding to the parameters of the chirp.
- For a pair of linear chirps, "par" should be {"SNR":[x_1, x_2], "a":[y_1, y_2], "b":z} where x_1 and x_2 are positive floats corresponding to the SNR of, respectively, the first and second chirp; y_1 and y_2 are floats corresponding to the intercept of the main axis of, respectively, the first and second chirp; and z is a float corresponding to the common slope of both chirps.

The optional arguments are:

- win_length: int -> Length of the window vector. By default it is 512.
- T0: float -> Time at which the signal starts. By default it is 0 for chirps and -T_length/2 for Hermite functions.
- T_length: float -> Time length of the signal. By default it is 16.
- fs: float -> Sampling frequency. By default it is 1.0.

The two main functions used for simulations are "Signal.plot" and "Signal.add_yellow_box"

Instructions for the function Signal.plot
~~~~~~

The full syntax of the function is the following

plot(show_max=False, plot_type=["spectro", "intensity"])

The main arguments are:

- show_max: bool -> If the local maxima should be drawn on the noisy spectrogram or not.
- plot_type: str list -> What type of plots should be drawn: "spectro" for the spectrogram of the noisy signal, "intensity" for the intensity of the zeros of the noisy spectrogram and "theo" for the spectrogram of the pure signal.

The function returns a figure and an Axes object of a matplotlib plot.

Instructions for the function Signal.add_yellow_box
~~~~~~

The full syntax of the function is the following

add_yellow_box(self, ax, N0)

The main arguments are:

- ax: Axes -> Axes object of a matplotlib plot.
- N0: int -> The number of the rectangle where a zero is trapped between two chirps.

The function returns nothing. It directly modify the Axes parameter to add a yellow rectangle corresponding to a rectangular trapping region for one of the zeros.

Dependencies
------------

This project was made with Python 3.10 and depends on the following packages:

-  numpy 2.2
-  scipy 1.15
-  matplotlib 3.10
-  cmocean 4.0

The following dependencies are optional, and only needed for the video generation.

-  tqdm 4.67
-  opencv 4.12

