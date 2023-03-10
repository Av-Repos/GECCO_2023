#include "tabu_search.h"
#include "RandomPerm.h"
#include "EvaluateQAP.h"
#include "globals.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

char * Instance;
int PermuSize, Function, Repetition, TabuSize, MaxEvals, Symmetry, mode;
int ** flow_matrix;
int ** zero_flow_matrix;
int ** dist_matrix;
int ** zero_dist_matrix;
long int m_totalpq;
long int ** m_pq;
long int **** sum;
long int *** sum_exp;
int *** its;
long double ** values;

//Reads a matrix from the specified instance.
void readMatrix(int **mat,  FILE *fd){
    int i,j,a;

    for(i=0;i<PermuSize;i++){
        for(j=0;j<PermuSize;j++){
            fscanf(fd,"%d\t",&a);
            mat[i][j] = (int) a;
        }
    }
}

//Executes the Tabu Search.
void evalRandomSample(){

	  strct_optimo optimo;

    initialize_vals();

	  optimo.opt_permu=malloc(PermuSize*sizeof(int));

    RandomPerm(optimo.opt_permu);
    optimo.opt_fitness=EvaluateQAP(0, optimo.opt_permu);

    tabu_search(Function,&optimo,TabuSize,MaxEvals);
}

//Main function.
int main (int argc, char *argv[]) {

	  int i,j,p,q,sym_dist,sym_flow;
    FILE *fd;

    /*
     1.- Instance file
     2.- Repetition number (random seed)
     3.- Tabu list size
     4.- Number of solution evaluations
     5.- Header (1 yes)
     */
    if(argc != 6){
      printf("Incorrect format.\n");
      exit(1);
    }

    Instance = argv[1];
    Repetition = atoi(argv[2]);
    Function = 0;
    TabuSize = atoi(argv[3]);
    MaxEvals = atoi(argv[4]);

    srand(Repetition);

    fd=fopen(Instance,"r");

    if(fd==NULL){
      printf("The specified file doesn't exist.\n");
      exit(2);
    }

    fscanf(fd,"%d%*[^\n]",&PermuSize);

    values = malloc(3*sizeof(long double*));
    values[0] = (long double*) malloc(3*sizeof(long double));
    values[1] = (long double*) malloc(3*sizeof(long double));
    values[2] = (long double*) malloc(3*sizeof(long double));

    dist_matrix = malloc(PermuSize*sizeof(int*));
    zero_dist_matrix = malloc(PermuSize*sizeof(int*));
    flow_matrix = malloc(PermuSize*sizeof(int*));
    zero_flow_matrix = malloc(PermuSize*sizeof(int*));
    for ( i=0; i<PermuSize; i++) {
		    dist_matrix[i] = (int *)malloc(PermuSize*sizeof(int));
        zero_dist_matrix[i] = (int *)malloc(PermuSize*sizeof(int));
        flow_matrix[i] = (int *)malloc(PermuSize*sizeof(int));
        zero_flow_matrix[i] = (int *)malloc(PermuSize*sizeof(int));
	  }

    //Reads distance and flow matrices.
    readMatrix(dist_matrix, fd);
    readMatrix(flow_matrix, fd);

    fclose(fd);

    for(i=0; i<PermuSize; i++){
      for(j=0; j<PermuSize; j++){
        if(i!=j){
          zero_dist_matrix[i][j] = dist_matrix[i][j];
          zero_flow_matrix[i][j] = flow_matrix[i][j];
        }else{
          zero_dist_matrix[i][j] = 0;
          zero_flow_matrix[i][j] = 0;
        }
      }
    }

    //CHECK SYMMETRY
    sym_dist = 1;

    for(i=0;i<PermuSize-1;i++){
      for(j=i+1;j<PermuSize;j++){
        if(dist_matrix[i][j]!=dist_matrix[j][i]) sym_dist = 0;
      }
    }

    sym_flow = 1;

    for(p=0;p<PermuSize-1;p++){
      for(q=p+1;q<PermuSize;q++){
        if(flow_matrix[p][q]!=flow_matrix[q][p]) sym_flow = 0;
      }
    }

    if(!sym_dist && !sym_flow)
      Symmetry = 0; //Asymmetric.
    else if(sym_dist && !sym_flow)
      Symmetry = 1; //Semi-symmetric.
    else if(!sym_dist && sym_flow)
      Symmetry = 2; //Semi-symmetric.
    else
      Symmetry = 3; //Symmetric.

    //CALCULATE VALUES FOR EvaluateQAP_all
    m_totalpq = 0;

    for (p=0;p<PermuSize;p++){
        for (q=0;q<p;q++)
            m_totalpq += zero_flow_matrix[p][q];
        for (q=p+1;q<PermuSize;q++)
            m_totalpq += zero_flow_matrix[p][q];
    }

    m_pq = malloc(PermuSize*sizeof(long int *));

    for(i=0;i<PermuSize;i++){
        m_pq[i] = (long int*) malloc(PermuSize*sizeof(long int));
        m_pq[i][i] = 0;
    }

    for(i=0;i<PermuSize-1;i++){
      for(j=i+1;j<PermuSize;j++){
        m_pq[i][j] = m_totalpq;

        p=i;
        for (q=0;q<i;q++){
            m_pq[i][j] -= zero_flow_matrix[p][q];
        }
        for (q=i+1;q<PermuSize;q++){
            m_pq[i][j] -= zero_flow_matrix[p][q];
        }

        //p==j denenan,
        p=j;
        for (q=0;q<j;q++){
            m_pq[i][j] -= zero_flow_matrix[p][q];
        }
        for (q=j+1;q<PermuSize;q++){
            m_pq[i][j] -= zero_flow_matrix[p][q];
        }

        //eta q-rekin berdin:

        //q==i denenan,
        q=i;
        for (p=0;p<i;p++){
            m_pq[i][j] -= zero_flow_matrix[p][q];
        }
        for (p=i+1;p<j;p++){		//--> p==j aurreko pausuan kendu dugu
            m_pq[i][j] -= zero_flow_matrix[p][q];
        }
        for (p=j+1;p<PermuSize;p++){
            m_pq[i][j] -= zero_flow_matrix[p][q];
        }

        //q==j denenan,
        q=j;
        for (p=0;p<i;p++){
            m_pq[i][j] -= zero_flow_matrix[p][q];
        }
        for (p=i+1;p<j;p++){		//--> p==j aurreko pausuan kendu dugu
            m_pq[i][j] -= zero_flow_matrix[p][q];
        }
        for (p=j+1;p<PermuSize;p++){
            m_pq[i][j] -= zero_flow_matrix[p][q];
        }
        m_pq[j][i] = m_pq[i][j];
      }
    }

    //CALCULATE VALUES FOR EvaluateQAP_change_all
    its = malloc(PermuSize*sizeof(int **));
    for(i=0;i<PermuSize;i++){
      its[i] = malloc(PermuSize*sizeof(int *));
      for(j=0;j<PermuSize;j++){
        if(i!=j)
          its[i][j] = (int*)malloc((PermuSize-2)*sizeof(int));
        else
          its[i][j] = (int*)malloc((PermuSize-1)*sizeof(int));
        q=0;
        for(p=0;p<PermuSize;p++){
          if(p!=i && p!=j){
            its[i][j][q]=p;
            q++;
          }
        }
      }
    }

    sum = malloc(PermuSize*sizeof(long int***));

    for ( i=0; i<PermuSize; i++) {
      sum[i] = malloc(PermuSize*sizeof(long int**));
      for(j=0;j<PermuSize;j++){
        sum[i][j] = malloc(PermuSize*sizeof(long int*));
        for(p=0;p<PermuSize;p++){
          sum[i][j][p] = (long int*) malloc(2*sizeof(long int));
          sum[i][j][p][0] = 0;
          sum[i][j][p][1] = 0;
          for(q=0;q<PermuSize;q++){
              if(q!=i && q!=j && q!=p){
                  sum[i][j][p][0] += zero_flow_matrix[q][i];
                  sum[i][j][p][0] -= zero_flow_matrix[q][j];
                  sum[i][j][p][1] += zero_flow_matrix[i][q];
                  sum[i][j][p][1] -= zero_flow_matrix[j][q];
              }
          }
        }
      }
    }

    sum_exp = malloc(PermuSize*sizeof(long int**));

    for ( i=0; i<PermuSize; i++) {
      sum_exp[i] = malloc(PermuSize*sizeof(long int*));
      for(j=0;j<PermuSize;j++){
        sum_exp[i][j] = (long int*) malloc(2*sizeof(long int));
        sum_exp[i][j][0] = 0;
        sum_exp[i][j][1] = 0;
        for(p=0;p<PermuSize;p++){
            if(p!=i && p!=j){
                sum_exp[i][j][0] += zero_flow_matrix[p][i];
                sum_exp[i][j][0] -= zero_flow_matrix[p][j];
                sum_exp[i][j][1] += zero_flow_matrix[i][p];
                sum_exp[i][j][1] -= zero_flow_matrix[j][p];
            }
        }
      }
    }

    //Writes output header.
    if(atoi(argv[5]) == 1)
      printf("Instance,Algorithm,Repetition,Iteration,f,f1,f2,f3,f1_1,f1_2,f1_3,f2_1,f2_2,f2_3,f3_1,f3_2,f3_3,Permutation\n");

    //Random execution of the Tabu Search.
    evalRandomSample();
}