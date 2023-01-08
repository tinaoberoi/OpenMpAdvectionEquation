# Project 1

## Serial Implementation

File `my_advection_program.c` contains the implementation of Advection Pseducode.

Commands to run the code:
`gcc my_advection_program.c -o my_advection_program -lm`
`./my_advection_program 400 10000 1.0 1.0e6 5.0e-7 2.85e-7`

The code will generate a guassian.txt (containing the initial values) and timestamp.txt (containing the values for selected n).

## Plotting

The plots are present in the `plots` folder for timestamps:
[100, 200, 500, 1000, 1500, 5000, 8000, 9000, 9500]

![Graph at timestamp 0](plots/image0.png)
![Graph at timestamp 100](plots/image100.png)
![Graph at timestamp 200](plots/image200.png)
![Graph at timestamp 500](plots/image500.png)
![Graph at timestamp 1000](plots/image1000.png)
![Graph at timestamp 1500](plots/image1500.png)
![Graph at timestamp 5000](plots/image5000.png)
![Graph at timestamp 8000](plots/image8000.png)
![Graph at timestamp 9000](plots/image9000.png)
![Graph at timestamp 9500](plots/image9500.png)

A gif is present as output.gif file created using the above mentioned plots.
![Output Animation](output.gif)
