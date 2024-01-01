/******************************************************************************
 * FILE: mpi_hello.c
 * DESCRIPTION:
 *   MPI tutorial example code: Simple hello world program
 * AUTHOR: Blaise Barney
 * LAST REVISED: 03/05/10
 ******************************************************************************/
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "graph.hpp"

#define MASTER 0
#define N 10

bool IsValidVertex(int x)
{
   return 1 <= x && x <= (N * N);
}

int main(int argc, char *argv[])
{
   int numtasks, taskid, len;
   char hostname[MPI_MAX_PROCESSOR_NAME];

   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
   MPI_Get_processor_name(hostname, &len);
   // printf ("Hello from task %d on %s!\n", taskid, hostname);
   // std::cout << "Hello from task " << taskid << " on " << hostname << std::endl;
   if (taskid == MASTER)
   {
      printf("MASTER: Number of MPI tasks is: %d\n", numtasks);
      Graph my_graph;

      for (int i = 1; i <= N * N; i++)
      {
         my_graph.AddVertex(i);
      }

      for (int i = 1; i <= N * N; i++)
      {

         if (IsValidVertex(i-N))
         {
            my_graph.AddEdge(i, i-N);
         }
         if (IsValidVertex(i+N))
         {
            my_graph.AddEdge(i, i+N);
         }
         if (i%N)
         {
            if (IsValidVertex(i+1))
            {
               my_graph.AddEdge(i, i+1);
            }
            
         }
         if ((i-1)%N)
         {
            if (IsValidVertex(i-1))
            {
               my_graph.AddEdge(i, i-1);
            }
            
         }
         
      }

      my_graph.Print();
      std::cout << "\n====BFS===\n";
      my_graph.InitSingleBFS(13);
      my_graph.RunBFSToStable();
      my_graph.PrintSingleBFS();
   }
   MPI_Finalize();
}
