#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <assert.h>

int main(int argc, char *argv[] ){
    int N, NT;
    double L, T, u, v;
    
    N = atoi(argv[1]);
    printf("%d\n", N);
    NT = atoi(argv[2]);
    printf("%d\n", NT);
    L = atof(argv[3]);
    printf("%f\n", L);
    T = atof(argv[4]);
    printf("%f\n", T);
    u = atof(argv[5]);
    printf("%.10f\n", u);
    v = atof(argv[6]);
    printf("%.10f\n", v);

    double** C_n1 = (double**)malloc(N*sizeof(double*));
    for(int i =0; i<N; i++){
        C_n1[i] = (double*)malloc(N*sizeof(double));
    }
    double** C_n2 = (double**)malloc(N*sizeof(double*));
    for(int i =0; i<N; i++){
        C_n2[i] = (double*)malloc(N*sizeof(double));
    }
    double dx = L/N;
    double dt = T/NT;
    double temp = dx/sqrt(2*(u*u + v*v));
    int i, j;
    assert(dt <= temp);
    for(i =0; i<N; i++){
        for(j = 0; j<N; j++){
            double valx = pow(i*dx - L/2, 2); 
            double valy = pow(j*dx - L/2, 2);
            double ran = (valx + valy)/(2*pow(L/4, 2));
            double temp = exp(-ran);
            C_n1[i][j] = temp;
        }
    }

    FILE *fptr;
    fptr = fopen("guassian.txt","w");
    if(fptr == NULL)
    {
        printf("Error!");   
        exit(1);             
    }
    for(i =0; i<N; i++){
        for(j =0; j<N; j++){
            //printf("%f, ", C_n1[i][j]);
            fprintf(fptr,"%lf, ", C_n1[i][j]);
        }
        //printf("\n");
        fprintf(fptr, "\n");
    }
    
    double up, down, right, left;
    for(int n =0; n< NT; n++){
        for(i =0; i<N; i++){
            for(j =0; j<N; j++){
                //printf("i = %d, j = %d, n = %d", i, j, n);
                up = C_n1[(((i-1)%N) + N)% N][j];
                down = C_n1[(((i+1)%N) + N)% N][j];
                left =  C_n1[i][(((j-1)%N) + N) %N];
                right = C_n1[i][((j+1%N) + N) %N];
                C_n2[i][j] = (up+down+left+right)/4 ;
                C_n2[i][j] -= (dt/(2*dx)) * (u*(down - up) + v*(right - left));
                //swap(C_n1, C_n2);
            }
        }
        for(int x =0; x<N; x++){
            for(int y =0; y<N; y++){
                C_n1[x][y] = C_n2[x][y];
            }
        }
        FILE *fptr2;
        if (n == 9500){
            fptr2 = fopen("timestamp.txt","w");
            if(fptr2 == NULL)
            {
                printf("Error!");   
                exit(1);             
            }
            for(int x =0; x<N; x++){
                for(int y =0; y<N; y++){
                    fprintf(fptr2,"%e, ", C_n1[x][y]);
                }
                fprintf(fptr2, "\n");
            }
        }   
    }
    return 0;
}