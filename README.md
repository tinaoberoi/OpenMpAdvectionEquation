# Project 1

## Milestone 1

## Serial Implementation

File `my_advection_program.c` contains the implementation of Advection Pseducode.

Commands to run the code:
`gcc my_advection_program.c -o my_advection_program -lm`
`./my_advection_program 400 10000 1.0 1.0e6 5.0e-7 2.85e-7`

The code will generate a guassian.txt (containing the initial values) and timestamp.txt (containing the values for selected n).

## Plotting

The plots are present in the `plots` folder for timestamps:
[100, 200, 500, 1000, 1500, 5000, 8000, 9000, 9500]

![Graph at timestamp 0](./scripts/plots/image0.png)
![Graph at timestamp 100](./scripts/plots/image1.png)
![Graph at timestamp 200](./scripts/plots/image2.png)
![Graph at timestamp 500](./scripts/plots/image3.png)
![Graph at timestamp 1000](./scripts/plots/image4.png)
![Graph at timestamp 1500](./scripts/plots/image5.png)
![Graph at timestamp 5000](./scripts/plots/image6.png)
![Graph at timestamp 8000](./scripts/plots/image7.png)
![Graph at timestamp 9000](./scripts/plots/image8.png)
![Graph at timestamp 9500](./scripts/plots/image9.png)

A gif is present as output.gif file created using the above mentioned plots.
![Output Animation](output.gif)


## Milestone 2

- Details of your processor and compiler choice for any data

- A demonstration that you get the same answer with the parallel and serial versions. Is bitwise
reproducibility expected?

- Your (best) grind rate (timesteps/s) using the parameters from Milestone 1 but on a problem with
N=10,000. This will be the standard value that we compare with each other in class.

The best performance was observed at 
- num_of_threads = 8 for outer loop
- num_of_threads = 8 for inner loop
- scheduled = static 
Time taken to execute the code for N = 10,000 varoed a lot on the linux server. But the best time taken was 578 sec.
According to this time the 
grid rate = N*N*NT/time
grid rate = 

- A plot of a strong scaling analysis using the following parameters:
∗ N = 3200
∗ NT = 400
∗ L = 1.0
∗ T = 1.0e3
∗ u = 5.0e-7
∗ v = 2.85e-7

- A plot of a second strong scaling study, using the same parameters as above but with N=200
For schedule = guided, N = 200, algo = first_order
![Graph at od weak scaling](./milestone2/images/200N_first_order_guided.png)

For schedule = lax, N = 200, algo = lax
![Graph at od weak scaling](./milestone2/images/200N_lax_guided.png)

For schedule = guided, N = 200, algo = second_order
![Graph at od weak scaling](./milestone2/images/200N_second_order_guided.png)

For schedule = static, N = 200, algo = lax
![Graph at od weak scaling](./milestone2/images/200N_lax_static.png)

For schedule = static, N = 200, algo = second_order
![Graph at od weak scaling](./milestone2/images/200N_second_order_static.png)

For schedule = static, N = 200, algo = first_order
![Graph at od weak scaling](./milestone2/images/200N_first_order_static.png)


- A plot of a weak scaling analysis using the following parameters:
![Graph at od weak scaling](./milestone2/images/weak_scale.png)


## Plots
