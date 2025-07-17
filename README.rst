Python code accompanying the paper `On two fundamental properties of the zeros of
spectrograms of noisy signals <??????????????>`_ with 
`Rémi Bardenet <https://rbardenet.github.io/>`_.

The file "Signal.py" contains the "Signal" class with the function needed to simulate the signals and associated spectrograms. An exemple of utilisation of these function
is given in the file "Figures in the paper.ipynb", recreating the figures apppearing in the associated paper. The file "Generating_animations_video.py" contains some
additional code to generate small videos showing how the zeros and local maxima of a noisy pair of paralle chirps behave when both chirps gets closer and closer from each other.
An example of such video is shown in the file "Anim.avi".

The Signal class
------------------------------------

Syntax
~~~~~~

The full syntax for defining an object of the Signal class is

Signal(par, win_length=512, T0=None, T_length=16, fs=1)

The main argument is "par" corresponding to a dictionnary of the parameters for the intended signal:

- For Hermite functions, "par" should be {"SNR":x, "k":y} where x is a positive float corresponding to the SNR of the singal and y is a positive int corresponding
to the parameter of the Hermire function.
- For a single linear chirp, "par" should be {"SNR":x, "a":y, "b":z} where x is a positive float corresponding to the SNR of the signal and y and z are floats corresponding
to the parameter of the chirp.
- For a pair of linear chirps, "par" should be {"SNR":[x_1, x_2], "a":[y_1, y_2], "b":z} where x_1 and x_2 are positive floats corresponding to the SNR of, respectively,
the first and second chirp; y_1 and y_2 are floats corresponding to the intercept of the main axis of respectively, the first and second chirp; and z is a float
corresponding to the common slope of both chirps.

The optional arguments are:

- win_length: int -> Length of the window vector.
- T0: float -> Time at which the signal starts. By default it is 0 for chirps and -T_length/2 for Hermite functions
- T_length: float -> Time length of the signal.
- fs: float -> Sampling frequency

The two main functions for elements of this class used for simulations are "Signal.plot" and "Signal.add_yellow_box"

Instructions for the function Signal.plot
-----------------------------------------

Syntax
~~~~~~

The full syntax of the function is the following

plot(show_max=False, plot_type=["spectro", "intensity"])

The main arguments are:

- show_max: bool -> The observed point pattern, an object of class "ppp".
- plot_type: str list -> The DPP family that is fitted to the data, it has to be either "Gauss", "Bessel" or "Cauchy".

The function returns a 2x2 matrix corresponding to the Fisher Information matrix for the parameters (rho, alpha) of classical stationnary DPP families.

Instructions for the function Signal.add_yellow_box
-----------------------------------------

Syntax
~~~~~~

The full syntax of the function is the following

add_yellow_box(self, ax, N0)

The main arguments are:

- ax: bool -> The observed point pattern, an object of class "ppp".
- N0: str list -> The DPP family that is fitted to the data, it has to be either "Gauss", "Bessel" or "Cauchy".

The function returns a 2x2 matrix corresponding to the Fisher Information matrix for the parameters (rho, alpha) of classical stationnary DPP families.

Dependencies
------------

This project depends on the following packages:

-  numpy
-  scipy
-  matplotlib
-  cmocean

The following dependencies are optional, and only needed for the video generation.

-  tqdm
-  opencv

