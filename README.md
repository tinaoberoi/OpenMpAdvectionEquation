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

The params are as follows:

N NT L T u v algorithm(lax/first_order/second_order) num_threads

To run serial code :

`make MyAdvectionSerial`

`./my_advection_program_serial 1131 400 1.0 1.0e3 5.0e-7 2.85e-7 lax 8`

To run parallel code:

`make MyAdvection`

`./my_advection_program_parallel 10000 20000 1.0 1.0e6 5.0e-7 2.85e-7 first_order 8`

- Details of your processor and compiler choice for any data

Used the linux server <br />

CPU(s): 16 <br />
Threads per core: 2 <br />
complier used : gcc <br />

- A demonstration that you get the same answer with the parallel and serial versions. Is bitwise
reproducibility expected?<br />

The output of both parallel and serial matches, output present in `milestone2/parallel_metric` and `milestone2/serial_metric`
for all three algos `lax, first_order, second_order`.
- Your (best) grind rate (timesteps/s) using the parameters from Milestone 1 but on a problem with
N=10,000. This will be the standard value that we compare with each other in class.<br />

The best performance was observed at <br />
num_of_threads = 8 for outer loop <br />
num_of_threads = 8 for inner loop <br />
scheduled = static <br />
Time taken to execute the code for N = 10,000 varoed a lot on the linux server. But the best time taken was 578 sec.<br />
According to this time the <br />
grid rate = N*N*NT/time <br />
grid rate = 3,527,336,860.670194 cells update per sec.<br />

## Plots

- A plot of a strong scaling analysis using the following parameters:
??? N = 3200
??? NT = 400
??? L = 1.0
??? T = 1.0e3
??? u = 5.0e-7
??? v = 2.85e-7

For schedule = guided, N = 3200, algo = first_order
![Graph at od weak scaling](./milestone2/images/3200N_first_order_guided.png)

For schedule = guided, N = 3200, algo = lax
![Graph at od weak scaling](./milestone2/images/3200N_lax_guided.png)

For schedule = guided, N = 3200, algo = second_order
![Graph at od weak scaling](./milestone2/images/3200N_second_order_guided.png)

For schedule = static, N = 3200, algo = lax
![Graph at od weak scaling](./milestone2/images/3200N_lax_static.png)


For schedule = static, N = 3200, algo = second_order
![Graph at od weak scaling](./milestone2/images/3200N_second_order_static.png)

For schedule = static, N = 3200, algo = first_order
![Graph at od weak scaling](./milestone2/images/3200N_first_order_static.png)


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

Since the server was very slow and ureliable I have generated the weak scales for the following:

N        gridpoints      n

200       40000          1

283       80000          2

400       160000         4

565       320000         8

800       640000         16

1131      1280000.       32

1600      2560000        64

![Graph at od weak scaling](./milestone2/images/weak_scale.png)


## Milestone 3

To run :

```
mpicc mpi_lax_nonuniform.c -fopenmp -o test2 -lm

```
### Description of parallel startergy

The master node (RANK = 0) broadcasts the params N, NT, u, v, L to all the other ranks. <br/>
Each rank executes its own chunk, and in the end sends its to master to create the final array.

The matrix generation (for loop) for invdividual ranks is parallelised using OpenMp. For this the optimum
time was coming at `num_threads = 4`.

### Demonstration of correctness
- GIF for lax order using MPI <br/>
<img src="./milestone3/lax.gif" width="500" height="500"/>

<br/><br/>
<hr>
- GIF for first order using MPI <br/>
<img src="./milestone3/first_order.gif" width="500" height="500"/>
<br/><br/>
<hr>
- GIF for second order using MPI <br/>
<img src="./milestone3/second_order.gif" width="500" height="500"/>
<br/><br/>
<hr>

#### Performance Analysis

- Strong Scaling

No of RANKS | 1 | 4 | 16 |
--- | --- | --- | --- 
Log(Ranks) (log ranks base 4) | 0 | 1 | 2 |
Time (sec) | 7128.5 | 1765 | 502.9 |
Strong Scaling (Time/7128.5) | 1.09 | 4.01 | 14.67 |

![Graph at of strong scaling](./milestone3/strong_scale.png)

- Weak Scaling

No of RANKS | 1 | 4 | 16 
--- | --- | --- | --- 
Log (log ranks base 4) | 0 | 1 | 2 |
N | 4000 | 8000 | 16000 | 
Weak Scaling (sec) | 1893.471 | 1634.109 | 1757.468 |


![Graph at of strong scaling](./milestone3/weak_scale.png)

No of RANKS | 1 | 4 | 16 |
--- | --- | --- | --- 
Grind Rate in Billions | 0.3 | 1.1 | 4.0 |


### Support of non uniform velocity

Here I tried different u and v relationships but only one made sense. But most of them were either moving very slow or movinf very fast.
The two relationships that worked here were:
- Moving very slow
The function I used here is `u = u + u/i` and `v = v+ v/j`

![Graph at of strong scaling](./milestone3/output_non_uniform.gif)

- The function I used here is `u = u + i*10e-10` and `v = v + j*10e-10`

![Graph at of strong scaling](./milestone3/output_non_uniform2.gif)


ACK : NISHCHAY KARLE 