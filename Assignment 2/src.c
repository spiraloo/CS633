#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include "mpi.h"

//  mpirun -np P -f hostfile ./halo Px N <num_time_steps> <seed> <stencil>
    // #pragma prutor-mpi-args: -np 8 -ppn 4
//   #pragma prutor-mpi-sysargs: 4 262144 10 7 5
//  #pragma prutor-mpi-sysargs: 4 16777216 10 7 9
//  #pragma prutor-mpi-sysargs: 4 4194304 10 7 5
//  #pragma prutor-mpi-sysargs: 4 4194304 10 7 9

int main( int argc, char *argv[])
{

    MPI_Status status;
    double sTime, eTime, time, max_value;

    MPI_Init (&argc, &argv);
    int myrank, P; 
    int position = 0;
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    // input argument test
    if(myrank == 0)
    {
        
        // printf("=====================\n");
        // printf("Px = %s\n", argv[1]);
        // printf("N = %s\n", argv[2]);
        // printf("num_time_steps = %s\n", argv[3]);
        // printf("seed = %s\n", argv[4]);
        // printf("stencil = %s\n", argv[5]);
   
    }
    int i, j, z;
    int Px = atoi(argv[1]);
    int Py = (int)P/Px;
    int D = (int)(sqrt(atoi(argv[2]))); // dimension of each process grid
    double **arr = (double **)malloc(D * sizeof(double *));
    double **new_arr = (double **)malloc(D * sizeof(double *));
    for (i=0; i<D; i++)
    {
        arr[i] = (double *)malloc(D * sizeof(double));
        new_arr[i] = (double *)malloc(D * sizeof(double));
    }
    int x = 2;
    int num_time_steps = atoi(argv[3]);
    // send and receive buffers
    
    double *rbuf = (double *)malloc(D*x*sizeof(double));
    double *sbuf = (double *)malloc(D*x*sizeof(double));


    // receive array
    double **rarr_top = (double**)malloc(x*sizeof(double*));
    double **rarr_bot = (double**)malloc(x*sizeof(double*));
    double **rarr_left = (double**)malloc(D*sizeof(double*));
    double **rarr_right = (double**)malloc(D*sizeof(double*));
    for(i=0; i<x; i++)
    {
        rarr_top[i] = (double*)malloc(D*sizeof(double));
        rarr_bot[i] = (double*)malloc(D*sizeof(double));
    }
    for(i=0; i<D; i++)
    {
        rarr_left[i] = (double*)malloc(x*sizeof(double));
        rarr_right[i] = (double*)malloc(x*sizeof(double));
    }
    int row = myrank / Px;
    int col = myrank % Px;
    
    int seed = atoi(argv[4]);
    srand(seed*(myrank+10));
    int subrank = 0;
    //make a subcommunicator for each row
    MPI_Comm row_comm;
    MPI_Comm_split(MPI_COMM_WORLD, row, myrank, &row_comm);
    MPI_Comm_rank(row_comm, &subrank);
    double ** gather_array_top = NULL;
    double ** gather_array_bot = NULL;
    double ** scatter_array_top = NULL;
    double ** scatter_array_bot = NULL;
    if (subrank == 0){
        gather_array_top = (double **)malloc(x * Px *  sizeof(double*));
        gather_array_bot = (double **)malloc(x * Px *  sizeof(double*));
        scatter_array_top = (double **)malloc(x * Px *  sizeof(double*));
        scatter_array_bot = (double **)malloc(x * Px *  sizeof(double*));
        for (i=0; i<x*Px; i++)
        {
            gather_array_top[i] = (double *)malloc(D * sizeof(double));
            gather_array_bot[i] = (double *)malloc(D * sizeof(double));
            scatter_array_top[i] = (double *)malloc(D * sizeof(double));
            scatter_array_bot[i] = (double *)malloc(D * sizeof(double));
        }
   
    }


    // z = 0;
    for(z = 0; z<2; z++)
    {
        
    for (i=0; i<D; i++)
    {
        for (j=0; j<D; j++)
        {
            arr[i][j] = abs(rand()+(i*rand()+j*myrank))/100;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(row_comm);
    
    sTime = MPI_Wtime();
    num_time_steps = atoi(argv[3]);
    while(num_time_steps--)
    {

        if(col % 2 == 0)
        {
            if(col != Px-1)
            {
                //pack using send buffer
                position = 0;
                for(i=0; i<D; i++)
                {
                    MPI_Pack(&arr[i][D-x], x, MPI_DOUBLE, sbuf, D*x*sizeof(double), &position, MPI_COMM_WORLD);
                }
                //send to right ( rank + 1 )
                MPI_Send(sbuf, position , MPI_PACKED, myrank+1, myrank + 1, MPI_COMM_WORLD); 

                //receive from right ( rank + 1 )
                MPI_Recv(rbuf, D*x, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD, &status);

                //unpack the received buffer into rarr_right
                position = 0;
                for(i=0; i<D; i++)
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
                for(i=0; i<D; i++)
                {
                    MPI_Unpack(rbuf, D*x*sizeof(double), &position, &rarr_left[i][0], x, MPI_DOUBLE, MPI_COMM_WORLD);
                }

                //pack using send buffer
                position = 0;
                for(i=0; i<D; i++)
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
                for(i=0; i<D; i++)
                {
                    MPI_Unpack(rbuf, D*x*sizeof(double), &position, &rarr_left[i][0], x, MPI_DOUBLE, MPI_COMM_WORLD);
                }

                //pack using send buffer
                position = 0;
                for(i=0; i<D; i++)
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
                for(i=0; i<D; i++)
                {
                    MPI_Pack(&arr[i][D-x], x, MPI_DOUBLE, sbuf, D*x*sizeof(double), &position, MPI_COMM_WORLD);
                }
                //send to right ( rank + 1 )
                MPI_Send(sbuf, position , MPI_PACKED, myrank+1, myrank + 1, MPI_COMM_WORLD); 

                //receive from right ( rank + 1 )
                MPI_Recv(rbuf, D*x, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD, &status);

                //unpack the received buffer into rarr_right
                position = 0;
                for(i=0; i<D; i++)
                {
                    MPI_Unpack(rbuf, D*x*sizeof(double), &position, &rarr_right[i][0], x, MPI_DOUBLE, MPI_COMM_WORLD);
                }
            }
        }

        //gather top and bottom
        if(z == 0)
        {  
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
        }
        else
        {

            if (subrank == 0){
            MPI_Gather(&arr[0][0],D*x,MPI_DOUBLE,&gather_array_top[0][0],D*x,MPI_DOUBLE,0,row_comm);     
            MPI_Gather(&(arr[D-x][0]),D*x,MPI_DOUBLE,&gather_array_bot[0][0],D*x,MPI_DOUBLE,0,row_comm);
            if (row % 2 == 0){
                if(row!=Py-1)
                {
                    MPI_Send(&gather_array_bot[0][0],x*Px*D,MPI_DOUBLE,myrank+Px,myrank,MPI_COMM_WORLD);
                    MPI_Recv(&scatter_array_bot[0][0],x*Px*D,MPI_DOUBLE,myrank+Px,myrank+Px,MPI_COMM_WORLD,&status);
                }
                if(row!=0)
                {
                    MPI_Recv(&scatter_array_top[0][0],x*Px*D,MPI_DOUBLE,myrank-Px,myrank-Px,MPI_COMM_WORLD,&status);
                    MPI_Send(&gather_array_top[0][0],x*Px*D,MPI_DOUBLE,myrank-Px,myrank,MPI_COMM_WORLD);
                }
            }
            else{
                if(row!=0)
                {
                    MPI_Recv(&scatter_array_top[0][0],x*Px*D,MPI_DOUBLE,myrank-Px,myrank-Px,MPI_COMM_WORLD,&status);
                    MPI_Send(&gather_array_top[0][0],x*Px*D,MPI_DOUBLE,myrank-Px,myrank,MPI_COMM_WORLD);
                }
                if(row!=Py-1)
                {
                    MPI_Send(&gather_array_bot[0][0],x*Px*D,MPI_DOUBLE,myrank+Px,myrank,MPI_COMM_WORLD);
                    MPI_Recv(&scatter_array_bot[0][0],x*Px*D,MPI_DOUBLE,myrank+Px,myrank+Px,MPI_COMM_WORLD,&status);
                }

            }
            
            //Scatter top and bottom to the recieve arrays top and botton of the intranode processes
            MPI_Scatter(&scatter_array_top[0][0],D*x,MPI_DOUBLE,&rarr_top[0][0],D*x,MPI_DOUBLE,0,row_comm);
            MPI_Scatter(&scatter_array_bot[0][0],D*x,MPI_DOUBLE,&rarr_bot[0][0],D*x,MPI_DOUBLE,0,row_comm);
            }
            else {
                MPI_Gather(&arr[0][0],D*x,MPI_DOUBLE,NULL,D*x,MPI_DOUBLE,0,row_comm);
                MPI_Gather(&(arr[D-x][0]),D*x,MPI_DOUBLE,NULL,D*x,MPI_DOUBLE,0,row_comm);
                MPI_Scatter(NULL,D*x,MPI_DOUBLE,&rarr_top[0][0],D*x,MPI_DOUBLE,0,row_comm);
                MPI_Scatter(NULL,D*x,MPI_DOUBLE,&rarr_bot[0][0],D*x,MPI_DOUBLE ,0,row_comm);
            }
        }


        for (i=0; i<D; i++)
        {
            for (j=0; j<D; j++)
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
                    new_arr[i][j] = temp/(count*1.0);
                }
                else
                {
                    new_arr[i][j] = (arr[i-1][j] + arr[i-2][j] + arr[i+1][j] +arr[i+2][j] + arr[i][j-1] + arr[i][j-2] + arr[i][j+1] + arr[i][j+2] + arr[i][j])/9.0;
                }
            }
        }
        // copy new_arr to arr
        for (i=0; i<D; i++)
        {
            for (j=0; j<D; j++)
            {
                arr[i][j] = new_arr[i][j];
            }
        }
        
        //add mpi barrier
        MPI_Barrier(MPI_COMM_WORLD);
    
    }
    eTime = MPI_Wtime();
    time = eTime - sTime;
    MPI_Reduce(&time, &max_value, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (myrank == 0){
    if (z == 0){
        printf("Time without leader = %lf\n",max_value);
    }
    else {
        printf("Time with leader = %lf\n",max_value);
    }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  
    }
   
    if (myrank == 0) {
        printf("Data=%lf\n", arr[0][0]);
    }
    
    if (subrank == 0){
        free(gather_array_top);
        free(gather_array_bot);
        free(scatter_array_top);
        free(scatter_array_bot);
    }
    

    for (i=0; i<D; i++)
    {
        free(arr[i]);
        free(new_arr[i]);
    }
    free(arr);
    free(new_arr);
    free(rbuf);
    free(sbuf);
    free(rarr_top);
    free(rarr_bot);
    MPI_Finalize();
    return EXIT_SUCCESS;
}