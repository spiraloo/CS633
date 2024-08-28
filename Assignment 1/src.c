#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include "mpi.h"


int main( int argc, char *argv[])
{

    MPI_Status status;
    double sTime, eTime, time, max_value;

    MPI_Init (&argc, &argv);
    int myrank; 
    int position = 0;
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);

   
    
    int Px = atoi(argv[1]);
    int Py = 3;
    int D = (int)(sqrt(atoi(argv[2])));
    double **arr = (double **)malloc(D * sizeof(double *));
    double **new_arr = (double **)malloc(D * sizeof(double *));
    for (int i=0; i<D; i++)
    {
        arr[i] = (double *)malloc(D * sizeof(double));
        new_arr[i] = (double *)malloc(D * sizeof(double));
    }
    int stencil = atoi(argv[5]);
    int x = (stencil == 5) ? 1 : 2;
    int num_time_steps = atoi(argv[3]);
    // send and receive buffers
    
    double *rbuf = (double *)malloc(D*x*sizeof(double));
    double *sbuf = (double *)malloc(D*x*sizeof(double));


    // receive array
    double **rarr_top = (double**)malloc(x*sizeof(double*));
    double **rarr_bot = (double**)malloc(x*sizeof(double*));
    double **rarr_left = (double**)malloc(D*sizeof(double*));
    double **rarr_right = (double**)malloc(D*sizeof(double*));
    for(int i=0; i<x; i++)
    {
        rarr_top[i] = (double*)malloc(D*sizeof(double));
        rarr_bot[i] = (double*)malloc(D*sizeof(double));
    }
    for(int i=0; i<D; i++)
    {
        rarr_left[i] = (double*)malloc(x*sizeof(double));
        rarr_right[i] = (double*)malloc(x*sizeof(double));
    }
    int row = myrank / Px;
    int col = myrank % Px;
    
    int seed = atoi(argv[4]);
    srand(seed*(myrank+10));

    // initialize data
    for (int i=0; i<D; i++)
    {
        for (int j=0; j<D; j++)
        {
            arr[i][j] = abs(rand()+(i*rand()+j*myrank))/100;
        }
    }
   
    

    sTime = MPI_Wtime();
    while(num_time_steps--)
    {

        if(col % 2 == 0)
        {
            if(col != Px-1)
            {
                //pack using send buffer
                position = 0;
                for(int i=0; i<D; i++)
                {
                    MPI_Pack(&arr[i][D-x], x, MPI_DOUBLE, sbuf, D*x*sizeof(double), &position, MPI_COMM_WORLD);
                }
                //send to right ( rank + 1 )
                MPI_Send(sbuf, position , MPI_PACKED, myrank+1, myrank + 1, MPI_COMM_WORLD); 

                //receive from right ( rank + 1 )
                MPI_Recv(rbuf, D*x, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD, &status);

                //unpack the received buffer into rarr_right
                position = 0;
                for(int i=0; i<D; i++)
                {
                    MPI_Unpack(rbuf, D*x*sizeof(double), &position, &rarr_right[i][0], x, MPI_DOUBLE, MPI_COMM_WORLD);
                }
            }
            if(col != 0)
            {
                // receive from left ( rank - 1 )
                MPI_Recv(rbuf, D*x, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD, &status);

                //unpack the received buffer into rarr_left
                position = 0;
                for(int i=0; i<D; i++)
                {
                    MPI_Unpack(rbuf, D*x*sizeof(double), &position, &rarr_left[i][0], x, MPI_DOUBLE, MPI_COMM_WORLD);
                }

                //pack using send buffer
                position = 0;
                for(int i=0; i<D; i++)
                {
                    MPI_Pack(&arr[i][0], x, MPI_DOUBLE, sbuf, D*x*sizeof(double), &position, MPI_COMM_WORLD);
                }

                //send to left ( rank - 1 )
                MPI_Send(sbuf, position , MPI_PACKED, myrank-1, myrank - 1, MPI_COMM_WORLD);
            }
        }

        else 
        {
            if(col != 0)
            {
                // receive from left ( rank - 1 )
                MPI_Recv(rbuf, D*x, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD, &status);

                //unpack the received buffer into rarr_left
                position = 0;
                for(int i=0; i<D; i++)
                {
                    MPI_Unpack(rbuf, D*x*sizeof(double), &position, &rarr_left[i][0], x, MPI_DOUBLE, MPI_COMM_WORLD);
                }

                //pack using send buffer
                position = 0;
                for(int i=0; i<D; i++)
                {
                    MPI_Pack(&arr[i][0], x, MPI_DOUBLE, sbuf, D*x*sizeof(double), &position, MPI_COMM_WORLD);
                }

                //send to left ( rank - 1 )
                MPI_Send(sbuf, position , MPI_PACKED, myrank-1, myrank - 1, MPI_COMM_WORLD);
            }

            if(col != Px-1)
            {
                //pack using send buffer
                position = 0;
                for(int i=0; i<D; i++)
                {
                    MPI_Pack(&arr[i][D-x], x, MPI_DOUBLE, sbuf, D*x*sizeof(double), &position, MPI_COMM_WORLD);
                }
                //send to right ( rank + 1 )
                MPI_Send(sbuf, position , MPI_PACKED, myrank+1, myrank + 1, MPI_COMM_WORLD); 

                //receive from right ( rank + 1 )
                MPI_Recv(rbuf, D*x, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD, &status);

                //unpack the received buffer into rarr_right
                position = 0;
                for(int i=0; i<D; i++)
                {
                    MPI_Unpack(rbuf, D*x*sizeof(double), &position, &rarr_right[i][0], x, MPI_DOUBLE, MPI_COMM_WORLD);
                }
            }
        }

        if(row % 2 == 0)
        {
            if(row != Py-1)
            {
                //pack using send buffer
                position = 0;
                MPI_Pack(&arr[D-x][0], D*x, MPI_DOUBLE, sbuf, D*x*sizeof(double), &position, MPI_COMM_WORLD);
                    
                //send to bottom ( rank + Px )
                MPI_Send(sbuf, position , MPI_PACKED, myrank+Px, myrank + Px, MPI_COMM_WORLD); 

                //receive from bottom ( rank + Px )
                MPI_Recv(rbuf, D*x, MPI_DOUBLE, myrank+Px, myrank, MPI_COMM_WORLD, &status);

                //unpack the received buffer into rarr_bot
                position = 0;
                MPI_Unpack(rbuf, D*x*sizeof(double), &position, &rarr_bot[0][0], D*x, MPI_DOUBLE, MPI_COMM_WORLD);
            }
            if(row != 0)
            {
                // receive from top ( rank - Px )
                MPI_Recv(rbuf, D*x, MPI_DOUBLE, myrank-Px, myrank, MPI_COMM_WORLD, &status);

                //unpack the received buffer into rarr_top
                position = 0;
                MPI_Unpack(rbuf, D*x*sizeof(double), &position, &rarr_top[0][0], D*x, MPI_DOUBLE, MPI_COMM_WORLD);
                
                //pack using send buffer
                position = 0;
                MPI_Pack(&arr[0][0], D*x, MPI_DOUBLE, sbuf, D*x*sizeof(double), &position, MPI_COMM_WORLD);

                //send to top ( rank - Px )
                MPI_Send(sbuf, position , MPI_PACKED, myrank-Px, myrank - Px, MPI_COMM_WORLD);
            }

        }
        else
        {
            if(row != 0)
            {
                // receive from top ( rank - Px )
                MPI_Recv(rbuf, D*x, MPI_DOUBLE, myrank-Px, myrank, MPI_COMM_WORLD, &status);

                //unpack the received buffer into rarr_top
                position = 0;
                MPI_Unpack(rbuf, D*x*sizeof(double), &position, &rarr_top[0][0], D*x, MPI_DOUBLE, MPI_COMM_WORLD);
                
                //pack using send buffer
                position = 0;
                MPI_Pack(&arr[0][0], D*x, MPI_DOUBLE, sbuf, D*x*sizeof(double), &position, MPI_COMM_WORLD);

                //send to top ( rank - Px )
                MPI_Send(sbuf, position , MPI_PACKED, myrank-Px, myrank - Px, MPI_COMM_WORLD);
            }

            if(row != Py-1)
            {
                //pack using send buffer
                position = 0;
                MPI_Pack(&arr[D-x][0], D*x, MPI_DOUBLE, sbuf, D*x*sizeof(double), &position, MPI_COMM_WORLD);
                    
                //send to bottom ( rank + Px )
                MPI_Send(sbuf, position , MPI_PACKED, myrank+Px, myrank + Px, MPI_COMM_WORLD); 

                //receive from bottom ( rank + Px )
                MPI_Recv(rbuf, D*x, MPI_DOUBLE, myrank+Px, myrank, MPI_COMM_WORLD, &status);

                //unpack the received buffer into rarr_bot
                position = 0;
                MPI_Unpack(rbuf, D*x*sizeof(double), &position, &rarr_bot[0][0], D*x, MPI_DOUBLE, MPI_COMM_WORLD);
            }
        }
        // Stencil computation

        if(stencil == 5)
        {
            for (int i=0; i<D; i++)
            {
                for (int j=0; j<D; j++)
                {
                    if(i == 0 || i == D-1 || j == 0 || j == D-1)
                    {
                        int count = 1;
                        double temp = arr[i][j];
                        // check bottom
                        if(i != D -1)
                        {
                            temp += arr[i+1][j];
                            count++;
                        }
                        else if(row != Py-1)
                        {
                            temp += rarr_bot[0][j];
                            count++;
                        }
                        // check top 
                        if(i != 0)
                        {
                            temp += arr[i-1][j];
                            count++;
                        }
                        else if(row != 0)
                        {
                            temp += rarr_top[0][j];
                            count++;
                        }
                        // check right
                        if(j != D-1)
                        {
                            temp += arr[i][j+1];
                            count++;
                        }
                        else if(col != Px-1)
                        {
                            temp += rarr_right[i][0];
                            count++;
                        }
                        // check left
                        if(j != 0)
                        {
                            temp += arr[i][j-1];
                            count++;
                        }
                        else if(col != 0)
                        {
                            temp += rarr_left[i][0];
                            count++;
                        }
                        new_arr[i][j] = temp/(1.0*count);
                    }
                    else
                    {
                        new_arr[i][j] = (arr[i-1][j] + arr[i+1][j] + arr[i][j-1] + arr[i][j+1] + arr[i][j])/5.0;
                    }
                }
            }
        }
        else
        {
            for (int i=0; i<D; i++)
            {
                for (int j=0; j<D; j++)
                {
                    if(i == 0 || i == 1 || i == D-1 || i == D-2 || j == 0  || j == 1|| j == D-1 || j == D-2)
                    {
                        int count = 1;
                        double temp = arr[i][j];

                        // check bottom
                        if( i < D -1)
                        {
                            temp += arr[i+1][j];
                            count++;

                            if(i!=D-2)
                            {
                                temp += arr[i+2][j];
                                count++;
                            }
                            else if(row != Py-1)
                            {
                                temp += rarr_bot[0][j];
                                count++;
                            }
                        }
                        else if(row != Py-1)
                        {
                            temp += rarr_bot[0][j]+rarr_bot[1][j];
                            count+=2;
                        }

                        // check top
                        if(i > 0)
                        {
                            temp += arr[i-1][j];
                            count++;

                            if(i!=1)
                            {
                                temp += arr[i-2][j];
                                count++;
                            }
                            else if(row != 0)
                            {
                                temp += rarr_top[1][j];
                                count+=1;
                            }
                        }
                        else if(row != 0)
                        {
                            temp += rarr_top[0][j]+rarr_top[1][j];
                            count+=2;
                        }

                        // check right
                        if(j < D-1)
                        {
                            temp += arr[i][j+1];
                            count++;

                            if(j!=D-2)
                            {
                                temp += arr[i][j+2];
                                count++;
                            }
                            else if(col != Px-1)
                            {
                                temp += rarr_right[i][0];
                                count+=1;
                            }
                        }
                        else if(col != Px-1)
                        {
                            temp += rarr_right[i][0]+rarr_right[i][1];
                            count+=2;
                        }

                        // check left
                        if(j > 0)
                        {
                            temp += arr[i][j-1];
                            count++;

                            if(j!=1)
                            {
                                temp += arr[i][j-2];
                                count++;
                            }
                            else if(col != 0)
                            {
                                temp += rarr_left[i][1];
                                count+=1;
                            }
                        }
                        else if(col != 0)
                        {
                            temp += rarr_left[i][0]+rarr_left[i][1];
                            count+=2;
                        }

                    }
                    else
                    {
                        new_arr[i][j] = (arr[i-1][j] + arr[i-2][j] + arr[i+1][j] +arr[i+2][j] + arr[i][j-1] + arr[i][j-2] + arr[i][j+1] + arr[i][j+2] + arr[i][j])/9.0;
                    }
                }
            }
        }
        // copy new_arr to arr
                for (int i=0; i<D; i++)
        {
            for (int j=0; j<D; j++)
            {
                arr[i][j] = new_arr[i][j];
            }
        }
      
        
        //mpi barrier
        MPI_Barrier(MPI_COMM_WORLD);
    
    }

    eTime = MPI_Wtime();
    time = eTime - sTime;

  
   
    MPI_Reduce(&time, &max_value, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (myrank == 0) {
        printf("Time Taken : %lf", max_value);
    }
    MPI_Finalize();
    free(arr);
    free(new_arr);
    free(rbuf);
    free(sbuf);
    free(rarr_top);
    free(rarr_bot);
    free(rarr_left);
    free(rarr_right);
    return 0;
}