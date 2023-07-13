#include "EvaluateQAP.h"
#include "globals.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

long int ** dif_ij;

/***********************************************************************************************************************************************************************/
//The following functions compute the f, f^1, f^2 and f^3 values of any given solution together with the corresponding sub-function values (f^m_chi, f^m_omega, f^m_tau).
//The sub-function values are stored in the "values" global variable.

//Computes f.
long double  evaluateQAP_function0_aux(int *genes){
    int FactA, i ,j;
    long double fitness = 0.0;
    for (i=0;i<PermuSize;i++){
        FactA = genes[i];
        for (j=0;j<PermuSize;j++){
            fitness= fitness + dist_matrix[i][j] * flow_matrix[FactA][genes[j]];
        }
    }
    return(fitness);
}

long double evaluateComp_function1_aux(int n, int i, int j, int p, int q, int genes_i, int genes_j){
  if(genes_i==p){
    if(genes_j==q){
       mode = 0;
       return n-3; //alfa
    }else{
       mode = 1;
       return -2;
    } //gamma
  }else if(genes_i==q){
    if(genes_j==p){
       mode = 0;
       return 1-n; //beta
    }else{
      mode = 1;
      return 0;
    } //delta
  }else if(genes_j==p){
    mode = 1;
    return 0; //delta
  }else if(genes_j==q){
    mode = 1;
    return -2; //gamma
  }else{
    mode = 2;
    return -1; //epsilon
  }
}

//Computes f^1, f^1_chi, f^1_omega and f^1_tau.
long double evaluateQAP_function1_aux(int *genes){
    int n, i, j, p, q;
    long double fitness, distAB,aux;

    n = PermuSize;
    fitness = 0.0;
    values[0][0] = 0.0;
    values[0][1] = 0.0;
    values[0][2] = 0.0;

	  for( i = 0; i < n;i++){
		    for( j = 0; j < n; j++){
             //if(i==j) continue;
             distAB = zero_dist_matrix[i][j];
             if(abs(distAB) < 0.000001) continue;
			       for( p = 0; p < n; p++){
				           for( q = 0; q < n; q++){
                        //if(p==q || flow_matrix[p][q] == 0) continue;
                        aux = evaluateComp_function1_aux(n, i, j, p, q, genes[i], genes[j]) * distAB * zero_flow_matrix[p][q];
                        values[0][mode] += aux;
                        fitness += aux;
				           }
             }
		    }
	  }

    values[0][0] /= (2*n);
    values[0][1] /= (2*n);
    values[0][2] /= (2*n);

    fitness = fitness/(2*n);
    return (fitness);
}

long double evaluateComp_function2_aux(int n, int i, int j, int p, int q, int genes_i, int genes_j){
  if(genes_i==p){
    if(genes_j==q){
       mode = 0;
       return n-3; //alfa
    }else{
       mode = 1;
       return 0;
    } //gamma
  }else if(genes_i==q){
    if(genes_j==p){
       mode = 0;
       return n-3; //beta
    }else{
      mode = 1;
      return 0;
    } //delta
  }else if(genes_j==p){
    mode = 1;
    return 0; //delta
  }else if(genes_j==q){
    mode = 1;
    return 0; //gamma
  }else{
    mode = 2;
    return 1; //epsilon
  }
}

//Computes f^2, f^2_chi, f^2_omega and f^2_tau.
long double evaluateQAP_function2_aux(int *genes){
    int n, i, j, p, q;
    long double fitness, distAB,aux;

    n = PermuSize;
    fitness = 0.0;
    values[1][0] = 0.0;
    values[1][1] = 0.0;
    values[1][2] = 0.0;

    for( i = 0; i < n;i++){
        for( j = 0; j < n; j++){
             //if(i==j) continue;
             distAB = zero_dist_matrix[i][j];
             if(abs(distAB) < 0.000001) continue;
             for( p = 0; p < n; p++){
                   for( q = 0; q < n; q++){
                        //if(p==q || flow_matrix[p][q] == 0) continue;
                        aux = evaluateComp_function2_aux(n, i, j, p, q, genes[i], genes[j]) * distAB * zero_flow_matrix[p][q];
                        values[1][mode] += aux;
                        fitness += aux;
                   }
             }
        }
    }

    values[1][0] /= (2*(n-2));
    values[1][1] /= (2*(n-2));
    values[1][2] /= (2*(n-2));

    fitness = fitness/(2*(n-2));
    return (fitness);
}

long double evaluateComp1_function3_aux(int n, int i, int j, int p, int q, int genes_i, int genes_j){
  if(genes_i==p){
    if(genes_j==q){
       mode = 0;
       return 2*n-3; //alfa
    }else{
       mode = 1;
       return n-2;
    } //gamma
  }else if(genes_i==q){
    if(genes_j==p){
       mode = 0;
       return 1; //beta
    }else{
      mode = 1;
      return 0;
    } //delta
  }else if(genes_j==p){
    mode = 1;
    return 0; //delta
  }else if(genes_j==q){
    mode = 1;
    return n-2; //gamma
  }else{
    mode = 2;
    return -1; //epsilon
  }
}

long double evaluateComp2_function3_aux(int i, int p, int genes_i){
	if(genes_i == p) return 1;
	else return 0;
}

//Computes f^3, f^3_chi, f^3_omega and f^3_tau.
long double evaluateQAP_function3_aux(int *genes){
    int n, i, j, p, q;
    long double fitness, distAB,aux;

    n = PermuSize;
    fitness = 0.0;
    values[2][0] = 0.0;
    values[2][1] = 0.0;
    values[2][2] = 0.0;

    for( i = 0; i < n;i++){
        for( j = 0; j < n; j++){
             //if(i==j) continue;
             distAB = zero_dist_matrix[i][j];
             if(abs(distAB) < 0.000001) continue;
             for( p = 0; p < n; p++){
                   for( q = 0; q < n; q++){
                        //if(p==q || flow_matrix[p][q] == 0) continue;
                        aux = evaluateComp1_function3_aux(n, i, j, p, q, genes[i], genes[j]) * distAB * zero_flow_matrix[p][q];
                        values[2][mode] += aux;
                        fitness += aux;
                   }
             }
        }
    }

    values[2][0] /= (n*(n-2));
    values[2][1] /= (n*(n-2));
    values[2][2] /= (n*(n-2));

    fitness = fitness/(n*(n-2));

    for (i=0;i<n;i++){
        distAB = dist_matrix[i][i];
        if(abs(distAB) < 0.000001) continue;
        for( p = 0; p < n; p++){
            //if(flow_matrix[p][p] == 0) continue;
            fitness = fitness + evaluateComp2_function3_aux(i, p, genes[i]) * distAB * flow_matrix[p][p];
        }
    }

    return (fitness);
}

//Auxiliary function that can be called outside the current file. It computes f, f^1, f^2 or f^3 (together with the corresponding sub-functions).
long double EvaluateQAP_aux(int fitnessFunctionId, int * genes){
    long double fitness = INT_MAX;
    switch (fitnessFunctionId){
        case 1:
            fitness = evaluateQAP_function1_aux(genes);
            break;
        case 2:
            fitness = evaluateQAP_function2_aux(genes);
            break;
        case 3:
            fitness = evaluateQAP_function3_aux(genes);
            break;
        default:
            fitness = evaluateQAP_function0_aux(genes);
            break;
    }

    return fitness;
}

/***********************************************************************************************************************************************************************/
//The following functions efficiently compute f, f^1, f^2 and f^3. The code has been optimized in order to allow a higher number of solution evaluations.

//Initializes the data structure that is needed to efficiently calculate the elementary function values.
void initialize_vals(){
  int i;
  dif_ij=(long int**)malloc(PermuSize*sizeof(long int*));
  for(i=0;i<PermuSize;i++){
    dif_ij[i]=(long int*)malloc(PermuSize*sizeof(long int));
  }
}

//Pre-computes the information that is needed to efficiently calculate the elementary function values.
void calculate_vals(int *genes){
    int i, j, k;
    for(i=0;i<PermuSize;i++){
      dif_ij[i][i] = 0;
      for(j=i+1;j<PermuSize;j++){
        dif_ij[i][j]=(dist_matrix[i][i]-dist_matrix[j][j])*(flow_matrix[genes[j]][genes[j]]-flow_matrix[genes[i]][genes[i]])+(dist_matrix[i][j]-dist_matrix[j][i])*(flow_matrix[genes[j]][genes[i]]-flow_matrix[genes[i]][genes[j]]);
        for(k=0;k<PermuSize;k++){
          if(k!=i && k!=j)
            dif_ij[i][j]+= (dist_matrix[i][k]-dist_matrix[j][k])*(flow_matrix[genes[j]][genes[k]]-flow_matrix[genes[i]][genes[k]])+(dist_matrix[k][i]-dist_matrix[k][j])*(flow_matrix[genes[k]][genes[j]]-flow_matrix[genes[k]][genes[i]]);
        }
      }
    }
}

//Re-computes the information that is needed to efficiently calculate the elementary function values.
void recalculate_vals(int *genes, int gene_1, int gene_2){
  int i, j, k;
  for(i=0;i<PermuSize;i++){
    for(j=i+1;j<PermuSize;j++){
      if(i!=gene_1&&i!=gene_2&&j!=gene_1&&j!=gene_2){
        dif_ij[i][j]+= (dist_matrix[gene_1][i]-dist_matrix[gene_1][j]+dist_matrix[gene_2][j]-dist_matrix[gene_2][i])*(flow_matrix[genes[gene_2]][genes[i]]-flow_matrix[genes[gene_2]][genes[j]]+flow_matrix[genes[gene_1]][genes[j]]-flow_matrix[genes[gene_1]][genes[i]]);
        dif_ij[i][j]+= (dist_matrix[i][gene_1]-dist_matrix[j][gene_1]+dist_matrix[j][gene_2]-dist_matrix[i][gene_2])*(flow_matrix[genes[i]][genes[gene_2]]-flow_matrix[genes[j]][genes[gene_2]]+flow_matrix[genes[j]][genes[gene_1]]-flow_matrix[genes[i]][genes[gene_1]]);
      }else{
        dif_ij[i][j]=(dist_matrix[i][i]-dist_matrix[j][j])*(flow_matrix[genes[j]][genes[j]]-flow_matrix[genes[i]][genes[i]])+(dist_matrix[i][j]-dist_matrix[j][i])*(flow_matrix[genes[j]][genes[i]]-flow_matrix[genes[i]][genes[j]]);
        for(k=0;k<PermuSize;k++){
          if(k!=i && k!=j)
            dif_ij[i][j]+= (dist_matrix[i][k]-dist_matrix[j][k])*(flow_matrix[genes[j]][genes[k]]-flow_matrix[genes[i]][genes[k]])+(dist_matrix[k][i]-dist_matrix[k][j])*(flow_matrix[genes[k]][genes[j]]-flow_matrix[genes[k]][genes[i]]);
        }
      }
    }
  }
}

//Computes f.
long double  evaluateQAP_function0(int *genes){
    int FactA, FactB, i ,j;
    long int fitness = 0;
    for (i=0;i<PermuSize;i++){
        FactA = genes[i];
        for (j=0;j<PermuSize;j++){
            FactB = genes[j];
            fitness= fitness + dist_matrix[i][j] * flow_matrix[FactA][FactB];
        }
    }
    return(fitness);
}

//Computes the f value difference between two neighbor solutions.
long double  evaluateQAP_change_function0(long double fitness, int gene_1, int gene_2){
    return(fitness+dif_ij[gene_1][gene_2]);
}


//Efficiently computes f^1.
long double evaluateQAP_function1_opt(int *genes){
    int i,j,p,q;
    int m_size = PermuSize;
    long int result=0;
    long int aux;
    int n3=m_size-3;
    int n1=1-m_size;
    int min,max;//,min_plus,max_plus;
    for (i=0;i<m_size;i++){
        for (j=0;j<i;j++){
            if(abs(zero_dist_matrix[i][j])<0.000001) continue;
            min=MIN(genes[i],genes[j]);
            max=MAX(genes[i],genes[j]);
            //min_plus=min+1;
            //max_plus=max+1;
            aux=-m_totalpq;

            //p==min denenan,
            p=min;
            for (q=0;q<min;q++)
                aux += zero_flow_matrix[p][q];
            for (q=min+1;q<m_size;q++)
                aux += zero_flow_matrix[p][q];

            //p==max denenan,
            p=max;
            for (q=0;q<max;q++)
                aux += zero_flow_matrix[p][q];
            for (q=max+1;q<m_size;q++)
                aux += zero_flow_matrix[p][q];

            //eta q-rekin berdin:

            //q==min denenan,
            q=min;
            for (p=0;p<min;p++)
                aux += zero_flow_matrix[p][q];
            for (p=min+1;p<max;p++)		//--> p==max aurreko pausuan kendu dugu
                aux += zero_flow_matrix[p][q];
            for (p=max+1;p<m_size;p++)
                aux += zero_flow_matrix[p][q];

            //q==max denenan,
            q=max;
            for (p=0;p<min;p++)
                aux += zero_flow_matrix[p][q];
            for (p=min+1;p<max;p++)		//--> p==max aurreko pausuan kendu dugu
                aux += zero_flow_matrix[p][q];
            for (p=max+1;p<m_size;p++)
                aux += zero_flow_matrix[p][q];

            for(q=0;q<genes[j];q++)
                aux += -2*zero_flow_matrix[genes[i]][q];
            for(q=genes[j]+1;q<m_size;q++)
                aux += -2*zero_flow_matrix[genes[i]][q];
            for(p=0;p<genes[i];p++)
                aux += -2*zero_flow_matrix[p][genes[j]];
            for(p=genes[i]+1;p<m_size;p++)
                aux += -2*zero_flow_matrix[p][genes[j]];

            //alpha case
            aux+=n3*zero_flow_matrix[genes[i]][genes[j]];

            //beta case.
            aux+=n1*zero_flow_matrix[genes[j]][genes[i]];

            //we multiply with distance here because the distance is a common factor in the calculation.
            result+=aux*zero_dist_matrix[i][j];
        }
        for (j=i+1;j<m_size;j++){
            if(abs(zero_dist_matrix[i][j])<0.000001) continue;
            min=MIN(genes[i],genes[j]);
            max=MAX(genes[i],genes[j]);
            //min_plus=min+1;
            //max_plus=max+1;
            aux=-m_totalpq;

            //p==min denenan,
            p=min;
            for (q=0;q<min;q++)
                aux += zero_flow_matrix[p][q];
            for (q=min+1;q<m_size;q++)
                aux += zero_flow_matrix[p][q];

            //p==max denenan,
            p=max;
            for (q=0;q<max;q++)
                aux += zero_flow_matrix[p][q];
            for (q=max+1;q<m_size;q++)
                aux += zero_flow_matrix[p][q];

            //eta q-rekin berdin:

            //q==min denenan,
            q=min;
            for (p=0;p<min;p++)
                aux += zero_flow_matrix[p][q];
            for (p=min+1;p<max;p++)		//--> p==max aurreko pausuan kendu dugu
                aux += zero_flow_matrix[p][q];
            for (p=max+1;p<m_size;p++)
                aux += zero_flow_matrix[p][q];

            //q==max denenan,
            q=max;
            for (p=0;p<min;p++)
                aux += zero_flow_matrix[p][q];
            for (p=min+1;p<max;p++)		//--> p==max aurreko pausuan kendu dugu
                aux += zero_flow_matrix[p][q];
            for (p=max+1;p<m_size;p++)
                aux += zero_flow_matrix[p][q];

            for(q=0;q<genes[j];q++)
                aux += -2*zero_flow_matrix[genes[i]][q];
            for(q=genes[j]+1;q<m_size;q++)
                aux += -2*zero_flow_matrix[genes[i]][q];
            for(p=0;p<genes[i];p++)
                aux += -2*zero_flow_matrix[p][genes[j]];
            for(p=genes[i]+1;p<m_size;p++)
                aux += -2*zero_flow_matrix[p][genes[j]];

            //alpha case
            aux+=n3*zero_flow_matrix[genes[i]][genes[j]];

            //beta case.
            aux+=n1*zero_flow_matrix[genes[j]][genes[i]];

            //we multiply with distance here because the distance is a common factor in the calculation.
            result+=aux*zero_dist_matrix[i][j];
        }
    }
    //result=result/(double)(2*(m_size-2));
    return (long double)result/(2*m_size);
}

//Computes the f^1 value difference between two neighbor solutions.
long double evaluateQAP_change_function1_opt(long double fitness, int *genes, int *genes_swap, int gene_1, int gene_2){
  int n, i, j, p, q,n1,n2,n24;
  int *it;
  int i_act, j_act, p_act, q_act, distAB;
  long int change,aux;

  n = PermuSize;
  n1 = 1-n;
  n2 = n-2;
  n24 = 2*n-4;
  change = 0;

  it=(int*)malloc(n2*sizeof(int));
  j=0;
  for(i=0;i<n;i++){
      if(i!=gene_1&&i!=gene_2){
        it[j]=i;
        j++;
      }
  }

  //Alpha, beta, gamma, delta, epsilon.
  for( i = 0; i < n2;i++){
     i_act = it[i];
     j_act = gene_1;
     distAB = zero_dist_matrix[i_act][j_act];
     if(abs(distAB) > 0.000001){
         aux = 0;
         for( p = 0; p < i; p++){
            p_act = it[p];
            q_act = gene_1;
            //From gamma to epsilon.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
            q_act = gene_2;
            //From epsilon to gamma.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         for( p = i+1; p < n2; p++){
            p_act = it[p];
            q_act = gene_1;
            //From gamma to epsilon.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
            q_act = gene_2;
            //From epsilon to gamma.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         p_act = i_act;
         q_act = gene_1;
         //From alpha to gamma.
         aux += n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         q_act = gene_2;
         //From gamma to alpha.
         aux -= n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         for( q = 0; q < i; q++){
            q_act = it[q];
            p_act = gene_1;
            //From delta to epsilon.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
            p_act = gene_2;
            //From epsilon to delta.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         for( q = i+1; q < n2; q++){
            q_act = it[q];
            p_act = gene_1;
            //From delta to epsilon.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
            p_act = gene_2;
            //From epsilon to delta.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         q_act = i_act;
         p_act = gene_1;
         //From beta to delta.
         aux -= n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         //From delta to beta.
         aux += n1*zero_flow_matrix[genes[p_act]][genes[q_act]];

         p_act = gene_1;
         q_act = gene_2;
         //From delta to gamma.
         aux -= 2*zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         q_act = gene_1;
         //From gamma to delta.
         aux += 2*zero_flow_matrix[genes[p_act]][genes[q_act]];
         change += aux*distAB;
     }
     j_act = gene_2;
     distAB = zero_dist_matrix[i_act][j_act];
     if(abs(distAB) > 0.000001){
         aux = 0;
         for( p = 0; p < i; p++){
            p_act = it[p];
            q_act = gene_1;
            //From epsilon to gamma.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
            q_act = gene_2;
            //From gamma to epsilon.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         for( p = i+1; p < n2; p++){
            p_act = it[p];
            q_act = gene_1;
            //From epsilon to gamma.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
            q_act = gene_2;
            //From gamma to epsilon.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         p_act = i_act;
         q_act = gene_1;
         //From gamma to alpha.
         aux -= n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         q_act = gene_2;
         //From alpha to gamma.
         aux += n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         for( q = 0; q < i; q++){
            q_act = it[q];
            p_act = gene_1;
            //From epsilon to delta.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
            p_act = gene_2;
            //From delta to epsilon.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         for( q = i+1; q < n2; q++){
            q_act = it[q];
            p_act = gene_1;
            //From epsilon to delta.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
            p_act = gene_2;
            //From delta to epsilon.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         q_act = i_act;
         p_act = gene_1;
         //From delta to beta.
         aux += n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         //From beta to delta.
         aux -= n1*zero_flow_matrix[genes[p_act]][genes[q_act]];

         p_act = gene_1;
         q_act = gene_2;
         //From gamma to delta.
         aux += 2*zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         q_act = gene_1;
         //From delta to gamma.
         aux -= 2*zero_flow_matrix[genes[p_act]][genes[q_act]];
         change += aux*distAB;
      }
  }

  for( j = 0; j < n2;j++){
     j_act = it[j];
     i_act = gene_1;
     distAB = zero_dist_matrix[i_act][j_act];
     if(abs(distAB) > 0.000001){
         aux = 0;
         for( p = 0; p < j; p++){
            p_act = it[p];
            q_act = gene_1;
            //From delta to epsilon.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
            q_act = gene_2;
            //From epsilon to delta.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         for( p = j+1; p < n2; p++){
            p_act = it[p];
            q_act = gene_1;
            //From delta to epsilon.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
            q_act = gene_2;
            //From epsilon to delta.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         p_act = j_act;
         q_act = gene_1;
         //From beta to delta.
         aux -= n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         q_act = gene_2;
         //From delta to beta.
         aux += n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         for( q = 0; q < j; q++){
            q_act = it[q];
            p_act = gene_1;
            //From gamma to epsilon.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
            p_act = gene_2;
            //From epsilon to gamma.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         for( q = j+1; q < n2; q++){
            q_act = it[q];
            p_act = gene_1;
            //From gamma to epsilon.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
            p_act = gene_2;
            //From epsilon to gamma.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         q_act = j_act;
         p_act = gene_1;
         //From alpha to gamma.
         aux += n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         //From gamma to alpha.
         aux -= n1*zero_flow_matrix[genes[p_act]][genes[q_act]];

         p_act = gene_1;
         q_act = gene_2;
         //From gamma to delta.
         aux += 2*zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         q_act = gene_1;
         //From delta to gamma.
         aux -= 2*zero_flow_matrix[genes[p_act]][genes[q_act]];
         change += aux*distAB;
     }
     i_act = gene_2;
     distAB = zero_dist_matrix[i_act][j_act];
     if(abs(distAB) > 0.000001){
         aux = 0;
         for( p = 0; p < j; p++){
            p_act = it[p];
            q_act = gene_1;
            //From epsilon to delta.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
            q_act = gene_2;
            //From delta to epsilon.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         for( p = j+1; p < n2; p++){
            p_act = it[p];
            q_act = gene_1;
            //From epsilon to delta.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
            q_act = gene_2;
            //From delta to epsilon.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         p_act = j_act;
         q_act = gene_1;
         //From delta to beta.
         aux += n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         q_act = gene_2;
         //From beta to delta.
         aux -= n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         for( q = 0; q < j; q++){
            q_act = it[q];
            p_act = gene_1;
            //From epsilon to gamma.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
            p_act = gene_2;
            //From gamma to epsilon.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         for( q = j+1; q < n2; q++){
            q_act = it[q];
            p_act = gene_1;
            //From epsilon to gamma.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
            p_act = gene_2;
            //From gamma to epsilon.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         q_act = j_act;
         p_act = gene_1;
         //From gamma to alpha.
         aux -= n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         //From alpha to gamma.
         aux += n1*zero_flow_matrix[genes[p_act]][genes[q_act]];

         p_act = gene_1;
         q_act = gene_2;
         //From delta to gamma.
         aux -= 2*zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         q_act = gene_1;
         //From gamma to delta.
         aux += 2*zero_flow_matrix[genes[p_act]][genes[q_act]];
         change += aux*distAB;
      }
  }

  i_act = gene_1;
  j_act = gene_2;

  distAB = zero_dist_matrix[i_act][j_act];
  if(abs(distAB) > 0.000001){
      aux = 0;
      for( p = 0; p < n2; p++){
         p_act = it[p];
         q_act = gene_1;
         //From delta to gamma.
         aux -= 2*zero_flow_matrix[genes[p_act]][genes[q_act]];
         q_act = gene_2;
         //From gamma to delta.
         aux += 2*zero_flow_matrix[genes[p_act]][genes[q_act]];
      }
      p_act = j_act;
      q_act = gene_1;
      //From beta to alpha.
      aux += n24*zero_flow_matrix[genes[p_act]][genes[q_act]];
      for( q = 0; q < n2; q++){
         q_act = it[q];
         p_act = gene_1;
         //From gamma to delta.
         aux += 2*zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         //From delta to gamma.
         aux -= 2*zero_flow_matrix[genes[p_act]][genes[q_act]];
      }
      q_act = j_act;
      p_act = gene_1;
      //From alpha to beta.
      aux -= n24*zero_flow_matrix[genes[p_act]][genes[q_act]];
      change += aux*distAB;
   }

   i_act = gene_2;
   j_act = gene_1;
   distAB = zero_dist_matrix[i_act][j_act];
   if(abs(distAB) > 0.000001){
       aux = 0;
       for( p = 0; p < n2; p++){
          p_act = it[p];
          q_act = gene_1;
          //From gamma to delta.
          aux += 2*zero_flow_matrix[genes[p_act]][genes[q_act]];
          q_act = gene_2;
          //From delta to gamma.
          aux -= 2*zero_flow_matrix[genes[p_act]][genes[q_act]];
       }
       p_act = j_act;
       q_act = gene_2;
       //From beta to alpha.
       aux += n24*zero_flow_matrix[genes[p_act]][genes[q_act]];
       for( q = 0; q < n2; q++){
          q_act = it[q];
          p_act = gene_1;
          //From delta to gamma.
          aux -= 2*zero_flow_matrix[genes[p_act]][genes[q_act]];
          p_act = gene_2;
          //From gamma to delta.
          aux += 2*zero_flow_matrix[genes[p_act]][genes[q_act]];
       }
       q_act = j_act;
       p_act = gene_2;
       //From alpha to beta.
       aux -= n24*zero_flow_matrix[genes[p_act]][genes[q_act]];
       change += aux*distAB;
    }

    fitness += (long double)change/(2*n);

    free(it);
    return (fitness);
}


//Efficiently computes f^2.
long double evaluateQAP_function2_opt(int *genes){
    int i,j,p,q;
    int m_size = PermuSize;
    long int result=0;
    long int aux;
    int n3=m_size-3;
    int min,max;//,min_plus,max_plus;
    for (i=0;i<m_size;i++){
        for (j=0;j<i;j++){
            if(abs(zero_dist_matrix[i][j])<0.000001) continue;
            min=MIN(genes[i],genes[j]);
            max=MAX(genes[i],genes[j]);
            //min_plus=min+1;
            //max_plus=max+1;
            aux=m_totalpq;

            //p==min denenan,
            p=min;
            for (q=0;q<min;q++)
                aux -= zero_flow_matrix[p][q];
            for (q=min+1;q<m_size;q++)
                aux -= zero_flow_matrix[p][q];

            //p==max denenan,
            p=max;
            for (q=0;q<max;q++)
                aux -= zero_flow_matrix[p][q];
            for (q=max+1;q<m_size;q++)
                aux -= zero_flow_matrix[p][q];

            //eta q-rekin berdin:

            //q==min denenan,
            q=min;
            for (p=0;p<min;p++)
                aux -= zero_flow_matrix[p][q];
            for (p=min+1;p<max;p++)		//--> p==max aurreko pausuan kendu dugu
                aux -= zero_flow_matrix[p][q];
            for (p=max+1;p<m_size;p++)
                aux -= zero_flow_matrix[p][q];

            //q==max denenan,
            q=max;
            for (p=0;p<min;p++)
                aux -= zero_flow_matrix[p][q];
            for (p=min+1;p<max;p++)		//--> p==max aurreko pausuan kendu dugu
                aux -= zero_flow_matrix[p][q];
            for (p=max+1;p<m_size;p++)
                aux -= zero_flow_matrix[p][q];

            //alpha case
            aux+=n3*zero_flow_matrix[genes[i]][genes[j]];

            //beta case.
            aux+=n3*zero_flow_matrix[genes[j]][genes[i]];

            //we multiply with distance here because the distance is a common factor in the calculation.
            result+=aux*zero_dist_matrix[i][j];
        }
        for (j=i+1;j<m_size;j++){
            if(abs(zero_dist_matrix[i][j])<0.000001) continue;
            min=MIN(genes[i],genes[j]);
            max=MAX(genes[i],genes[j]);
            //min_plus=min+1;
            //max_plus=max+1;
            aux=m_totalpq;

            //p==min denenan,
            p=min;
            for (q=0;q<min;q++)
                aux -= zero_flow_matrix[p][q];
            for (q=min+1;q<m_size;q++)
                aux -= zero_flow_matrix[p][q];

            //p==max denenan,
            p=max;
            for (q=0;q<max;q++)
                aux -= zero_flow_matrix[p][q];
            for (q=max+1;q<m_size;q++)
                aux -= zero_flow_matrix[p][q];

            //eta q-rekin berdin:

            //q==min denenan,
            q=min;
            for (p=0;p<min;p++)
                aux -= zero_flow_matrix[p][q];
            for (p=min+1;p<max;p++)		//--> p==max aurreko pausuan kendu dugu
                aux -= zero_flow_matrix[p][q];
            for (p=max+1;p<m_size;p++)
                aux -= zero_flow_matrix[p][q];

            //q==max denenan,
            q=max;
            for (p=0;p<min;p++)
                aux -= zero_flow_matrix[p][q];
            for (p=min+1;p<max;p++)		//--> p==max aurreko pausuan kendu dugu
                aux -= zero_flow_matrix[p][q];
            for (p=max+1;p<m_size;p++)
                aux -= zero_flow_matrix[p][q];

            //alpha case
            aux+=n3*zero_flow_matrix[genes[i]][genes[j]];

            //beta case.
            aux+=n3*zero_flow_matrix[genes[j]][genes[i]];

            //we multiply with distance here because the distance is a common factor in the calculation.
            result+=aux*zero_dist_matrix[i][j];
        }
    }
    //result=result/(double)(2*(m_size-2));
    return (long double)result/(2*(m_size-2));
}

//Computes the f^2 value difference between two neighbor solutions.
long double evaluateQAP_change_function2_opt(long double fitness, int *genes, int *genes_swap, int gene_1, int gene_2){
  int n, i, j, p, q, n2,n3;
  int *it;
  int i_act, j_act, p_act, q_act, distAB;
  long int change,aux;

  n = PermuSize;
  n2 = n-2;
  n3 = n-3;
  change = 0;

  it=(int*)malloc(n2*sizeof(int));
  j=0;
  for(i=0;i<n;i++){
      if(i!=gene_1&&i!=gene_2){
        it[j]=i;
        j++;
      }
  }

  //Alpha, beta, gamma, delta, epsilon.
  for( i = 0; i < n2;i++){
     i_act = it[i];
     j_act = gene_1;
     distAB = zero_dist_matrix[i_act][j_act];
     if(abs(distAB) > 0.000001){
         aux = 0;
         for( p = 0; p < i; p++){
            p_act = it[p];
            q_act = gene_1;
            //From gamma to epsilon.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
            q_act = gene_2;
            //From epsilon to gamma.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         for( p = i+1; p < n2; p++){
            p_act = it[p];
            q_act = gene_1;
            //From gamma to epsilon.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
            q_act = gene_2;
            //From epsilon to gamma.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         p_act = i_act;
         q_act = gene_1;
         //From alpha to gamma.
         aux -= n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         q_act = gene_2;
         //From gamma to alpha.
         aux += n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         for( q = 0; q < i; q++){
            q_act = it[q];
            p_act = gene_1;
            //From delta to epsilon.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
            p_act = gene_2;
            //From epsilon to delta.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         for( q = i+1; q < n2; q++){
            q_act = it[q];
            p_act = gene_1;
            //From delta to epsilon.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
            p_act = gene_2;
            //From epsilon to delta.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         q_act = i_act;
         p_act = gene_1;
         //From beta to delta.
         aux -= n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         //From delta to beta.
         aux += n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         change += aux*distAB;
     }
     j_act = gene_2;
     distAB = zero_dist_matrix[i_act][j_act];
     if(abs(distAB) > 0.000001){
         aux = 0;
         for( p = 0; p < i; p++){
            p_act = it[p];
            q_act = gene_1;
            //From epsilon to gamma.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
            q_act = gene_2;
            //From gamma to epsilon.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         for( p = i+1; p < n2; p++){
            p_act = it[p];
            q_act = gene_1;
            //From epsilon to gamma.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
            q_act = gene_2;
            //From gamma to epsilon.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         p_act = i_act;
         q_act = gene_1;
         //From gamma to alpha.
         aux += n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         q_act = gene_2;
         //From alpha to gamma.
         aux -= n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         for( q = 0; q < i; q++){
            q_act = it[q];
            p_act = gene_1;
            //From epsilon to delta.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
            p_act = gene_2;
            //From delta to epsilon.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         for( q = i+1; q < n2; q++){
            q_act = it[q];
            p_act = gene_1;
            //From epsilon to delta.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
            p_act = gene_2;
            //From delta to epsilon.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         q_act = i_act;
         p_act = gene_1;
         //From delta to beta.
         aux += n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         //From beta to delta.
         aux -= n3*zero_flow_matrix[genes[p_act]][genes[q_act]];

         change += aux*distAB;
      }
  }

  for( j = 0; j < n2;j++){
     j_act = it[j];
     i_act = gene_1;
     distAB = zero_dist_matrix[i_act][j_act];
     if(abs(distAB) > 0.000001){
         aux = 0;
         for( p = 0; p < j; p++){
            p_act = it[p];
            q_act = gene_1;
            //From delta to epsilon.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
            q_act = gene_2;
            //From epsilon to delta.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         for( p = j+1; p < n2; p++){
            p_act = it[p];
            q_act = gene_1;
            //From delta to epsilon.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
            q_act = gene_2;
            //From epsilon to delta.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         p_act = j_act;
         q_act = gene_1;
         //From beta to delta.
         aux -= n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         q_act = gene_2;
         //From delta to beta.
         aux += n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         for( q = 0; q < j; q++){
            q_act = it[q];
            p_act = gene_1;
            //From gamma to epsilon.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
            p_act = gene_2;
            //From epsilon to gamma.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         for( q = j+1; q < n2; q++){
            q_act = it[q];
            p_act = gene_1;
            //From gamma to epsilon.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
            p_act = gene_2;
            //From epsilon to gamma.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         q_act = j_act;
         p_act = gene_1;
         //From alpha to gamma.
         aux -= n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         //From gamma to alpha.
         aux += n3*zero_flow_matrix[genes[p_act]][genes[q_act]];

         change += aux*distAB;
     }
     i_act = gene_2;
     distAB = zero_dist_matrix[i_act][j_act];
     if(abs(distAB) > 0.000001){
         aux = 0;
         for( p = 0; p < j; p++){
            p_act = it[p];
            q_act = gene_1;
            //From epsilon to delta.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
            q_act = gene_2;
            //From delta to epsilon.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         for( p = j+1; p < n2; p++){
            p_act = it[p];
            q_act = gene_1;
            //From epsilon to delta.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
            q_act = gene_2;
            //From delta to epsilon.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         p_act = j_act;
         q_act = gene_1;
         //From delta to beta.
         aux += n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         q_act = gene_2;
         //From beta to delta.
         aux -= n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         for( q = 0; q < j; q++){
            q_act = it[q];
            p_act = gene_1;
            //From epsilon to gamma.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
            p_act = gene_2;
            //From gamma to epsilon.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         for( q = j+1; q < n2; q++){
            q_act = it[q];
            p_act = gene_1;
            //From epsilon to gamma.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
            p_act = gene_2;
            //From gamma to epsilon.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         q_act = j_act;
         p_act = gene_1;
         //From gamma to alpha.
         aux += n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         //From alpha to gamma.
         aux -= n3*zero_flow_matrix[genes[p_act]][genes[q_act]];

         change += aux*distAB;
      }
  }

  fitness += (long double)change/(2*n2);

  free(it);
  return (fitness);
}

//Efficiently computes f^3.
long double evaluateQAP_function3_opt(int *genes){
    int i,j,p,q;
    int m_size = PermuSize;
    long int result1=0;
    long int result2=0;
    long int aux;
    int n23=2*m_size-3;
    int n2=m_size-2;
    int min,max;//,min_plus,max_plus;
    for (i=0;i<m_size;i++){
        for (j=0;j<i;j++){
            if(abs(zero_dist_matrix[i][j])<0.000001) continue;
            min=MIN(genes[i],genes[j]);
            max=MAX(genes[i],genes[j]);
            //min_plus=min+1;
            //max_plus=max+1;
            aux=-m_totalpq;

            //p==min denenan,
            p=min;
            for (q=0;q<min;q++)
                aux += zero_flow_matrix[p][q];
            for (q=min+1;q<m_size;q++)
                aux += zero_flow_matrix[p][q];

            //p==max denenan,
            p=max;
            for (q=0;q<max;q++)
                aux += zero_flow_matrix[p][q];
            for (q=max+1;q<m_size;q++)
                aux += zero_flow_matrix[p][q];

            //eta q-rekin berdin:

            //q==min denenan,
            q=min;
            for (p=0;p<min;p++)
                aux += zero_flow_matrix[p][q];
            for (p=min+1;p<max;p++)		//--> p==max aurreko pausuan kendu dugu
                aux += zero_flow_matrix[p][q];
            for (p=max+1;p<m_size;p++)
                aux += zero_flow_matrix[p][q];

            //q==max denenan,
            q=max;
            for (p=0;p<min;p++)
                aux += zero_flow_matrix[p][q];
            for (p=min+1;p<max;p++)		//--> p==max aurreko pausuan kendu dugu
                aux += zero_flow_matrix[p][q];
            for (p=max+1;p<m_size;p++)
                aux += zero_flow_matrix[p][q];

            for(q=0;q<genes[j];q++)
                aux += n2*zero_flow_matrix[genes[i]][q];
            for(q=genes[j]+1;q<m_size;q++)
                aux += n2*zero_flow_matrix[genes[i]][q];
            for(p=0;p<genes[i];p++)
                aux += n2*zero_flow_matrix[p][genes[j]];
            for(p=genes[i]+1;p<m_size;p++)
                aux += n2*zero_flow_matrix[p][genes[j]];

            //alpha case
            aux+=n23*zero_flow_matrix[genes[i]][genes[j]];

            //beta case.
            aux+=zero_flow_matrix[genes[j]][genes[i]];

            //we multiply with distance here because the distance is a common factor in the calculation.
            result1+=aux*zero_dist_matrix[i][j];
        }
        for (j=i+1;j<m_size;j++){
            if(abs(zero_dist_matrix[i][j])<0.000001) continue;
            min=MIN(genes[i],genes[j]);
            max=MAX(genes[i],genes[j]);
            //min_plus=min+1;
            //max_plus=max+1;
            aux=-m_totalpq;

            //p==min denenan,
            p=min;
            for (q=0;q<min;q++)
                aux += zero_flow_matrix[p][q];
            for (q=min+1;q<m_size;q++)
                aux += zero_flow_matrix[p][q];

            //p==max denenan,
            p=max;
            for (q=0;q<max;q++)
                aux += zero_flow_matrix[p][q];
            for (q=max+1;q<m_size;q++)
                aux += zero_flow_matrix[p][q];

            //eta q-rekin berdin:

            //q==min denenan,
            q=min;
            for (p=0;p<min;p++)
                aux += zero_flow_matrix[p][q];
            for (p=min+1;p<max;p++)		//--> p==max aurreko pausuan kendu dugu
                aux += zero_flow_matrix[p][q];
            for (p=max+1;p<m_size;p++)
                aux += zero_flow_matrix[p][q];

            //q==max denenan,
            q=max;
            for (p=0;p<min;p++)
                aux += zero_flow_matrix[p][q];
            for (p=min+1;p<max;p++)		//--> p==max aurreko pausuan kendu dugu
                aux += zero_flow_matrix[p][q];
            for (p=max+1;p<m_size;p++)
                aux += zero_flow_matrix[p][q];

            for(q=0;q<genes[j];q++)
                aux += n2*zero_flow_matrix[genes[i]][q];
            for(q=genes[j]+1;q<m_size;q++)
                aux += n2*zero_flow_matrix[genes[i]][q];
            for(p=0;p<genes[i];p++)
                aux += n2*zero_flow_matrix[p][genes[j]];
            for(p=genes[i]+1;p<m_size;p++)
                aux += n2*zero_flow_matrix[p][genes[j]];

            //alpha case
            aux+=n23*zero_flow_matrix[genes[i]][genes[j]];

            //beta case.
            aux+=zero_flow_matrix[genes[j]][genes[i]];

            //we multiply with distance here because the distance is a common factor in the calculation.
            result1+=aux*zero_dist_matrix[i][j];
        }
    }

    for (i=0;i<m_size;i++){
        result2 += dist_matrix[i][i]*flow_matrix[genes[i]][genes[i]];
    }

    return (long double)result1/(m_size*(m_size-2))+(long double)result2;
}

//Computes the f^3 value difference between two neighbor solutions.
long double evaluateQAP_change_function3_opt(long double fitness, int *genes, int *genes_swap, int gene_1, int gene_2){
  int n, i, j, p, q,n1,n2,n24;
  int *it;
  int i_act, j_act, p_act, q_act, distAB;
  long int change,aux;

  n = PermuSize;
  n1 = 1-n;
  n2 = n-2;
  n24 = 2*n-4;
  change = 0;

  it=(int*)malloc(n2*sizeof(int));
  j=0;
  for(i=0;i<n;i++){
      if(i!=gene_1&&i!=gene_2){
        it[j]=i;
        j++;
      }
  }

  //Alpha, beta, gamma, delta, epsilon.
  for( i = 0; i < n2;i++){
     i_act = it[i];
     j_act = gene_1;
     distAB = zero_dist_matrix[i_act][j_act];
     if(abs(distAB) > 0.000001){
         aux = 0;
         for( p = 0; p < i; p++){
            p_act = it[p];
            q_act = gene_1;
            //From gamma to epsilon.
            aux += n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
            q_act = gene_2;
            //From epsilon to gamma.
            aux -= n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         for( p = i+1; p < n2; p++){
            p_act = it[p];
            q_act = gene_1;
            //From gamma to epsilon.
            aux += n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
            q_act = gene_2;
            //From epsilon to gamma.
            aux -= n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         p_act = i_act;
         q_act = gene_1;
         //From alpha to gamma.
         aux += n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         q_act = gene_2;
         //From gamma to alpha.
         aux -= n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         for( q = 0; q < i; q++){
            q_act = it[q];
            p_act = gene_1;
            //From delta to epsilon.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
            p_act = gene_2;
            //From epsilon to delta.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         for( q = i+1; q < n2; q++){
            q_act = it[q];
            p_act = gene_1;
            //From delta to epsilon.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
            p_act = gene_2;
            //From epsilon to delta.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         q_act = i_act;
         p_act = gene_1;
         //From beta to delta.
         aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         //From delta to beta.
         aux += zero_flow_matrix[genes[p_act]][genes[q_act]];

         p_act = gene_1;
         q_act = gene_2;
         //From delta to gamma.
         aux += n2*zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         q_act = gene_1;
         //From gamma to delta.
         aux -= n2*zero_flow_matrix[genes[p_act]][genes[q_act]];
         change += aux*distAB;
     }
     j_act = gene_2;
     distAB = zero_dist_matrix[i_act][j_act];
     if(abs(distAB) > 0.000001){
         aux = 0;
         for( p = 0; p < i; p++){
            p_act = it[p];
            q_act = gene_1;
            //From epsilon to gamma.
            aux -= n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
            q_act = gene_2;
            //From gamma to epsilon.
            aux += n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         for( p = i+1; p < n2; p++){
            p_act = it[p];
            q_act = gene_1;
            //From epsilon to gamma.
            aux -= n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
            q_act = gene_2;
            //From gamma to epsilon.
            aux += n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         p_act = i_act;
         q_act = gene_1;
         //From gamma to alpha.
         aux -= n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         q_act = gene_2;
         //From alpha to gamma.
         aux += n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         for( q = 0; q < i; q++){
            q_act = it[q];
            p_act = gene_1;
            //From epsilon to delta.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
            p_act = gene_2;
            //From delta to epsilon.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         for( q = i+1; q < n2; q++){
            q_act = it[q];
            p_act = gene_1;
            //From epsilon to delta.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
            p_act = gene_2;
            //From delta to epsilon.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         q_act = i_act;
         p_act = gene_1;
         //From delta to beta.
         aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         //From beta to delta.
         aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];

         p_act = gene_1;
         q_act = gene_2;
         //From gamma to delta.
         aux -= n2*zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         q_act = gene_1;
         //From delta to gamma.
         aux += n2*zero_flow_matrix[genes[p_act]][genes[q_act]];
         change += aux*distAB;
      }
  }

  for( j = 0; j < n2;j++){
     j_act = it[j];
     i_act = gene_1;
     distAB = zero_dist_matrix[i_act][j_act];
     if(abs(distAB) > 0.000001){
         aux = 0;
         for( p = 0; p < j; p++){
            p_act = it[p];
            q_act = gene_1;
            //From delta to epsilon.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
            q_act = gene_2;
            //From epsilon to delta.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         for( p = j+1; p < n2; p++){
            p_act = it[p];
            q_act = gene_1;
            //From delta to epsilon.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
            q_act = gene_2;
            //From epsilon to delta.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         p_act = j_act;
         q_act = gene_1;
         //From beta to delta.
         aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
         q_act = gene_2;
         //From delta to beta.
         aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
         for( q = 0; q < j; q++){
            q_act = it[q];
            p_act = gene_1;
            //From gamma to epsilon.
            aux += n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
            p_act = gene_2;
            //From epsilon to gamma.
            aux -= n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         for( q = j+1; q < n2; q++){
            q_act = it[q];
            p_act = gene_1;
            //From gamma to epsilon.
            aux += n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
            p_act = gene_2;
            //From epsilon to gamma.
            aux -= n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         q_act = j_act;
         p_act = gene_1;
         //From alpha to gamma.
         aux += n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         //From gamma to alpha.
         aux -= n1*zero_flow_matrix[genes[p_act]][genes[q_act]];

         p_act = gene_1;
         q_act = gene_2;
         //From gamma to delta.
         aux -= n2*zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         q_act = gene_1;
         //From delta to gamma.
         aux += n2*zero_flow_matrix[genes[p_act]][genes[q_act]];
         change += aux*distAB;
     }
     i_act = gene_2;
     distAB = zero_dist_matrix[i_act][j_act];
     if(abs(distAB) > 0.000001){
         aux = 0;
         for( p = 0; p < j; p++){
            p_act = it[p];
            q_act = gene_1;
            //From epsilon to delta.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
            q_act = gene_2;
            //From delta to epsilon.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         for( p = j+1; p < n2; p++){
            p_act = it[p];
            q_act = gene_1;
            //From epsilon to delta.
            aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
            q_act = gene_2;
            //From delta to epsilon.
            aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         p_act = j_act;
         q_act = gene_1;
         //From delta to beta.
         aux += zero_flow_matrix[genes[p_act]][genes[q_act]];
         q_act = gene_2;
         //From beta to delta.
         aux -= zero_flow_matrix[genes[p_act]][genes[q_act]];
         for( q = 0; q < j; q++){
            q_act = it[q];
            p_act = gene_1;
            //From epsilon to gamma.
            aux -= n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
            p_act = gene_2;
            //From gamma to epsilon.
            aux += n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         for( q = j+1; q < n2; q++){
            q_act = it[q];
            p_act = gene_1;
            //From epsilon to gamma.
            aux -= n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
            p_act = gene_2;
            //From gamma to epsilon.
            aux += n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         }
         q_act = j_act;
         p_act = gene_1;
         //From gamma to alpha.
         aux -= n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         //From alpha to gamma.
         aux += n1*zero_flow_matrix[genes[p_act]][genes[q_act]];

         p_act = gene_1;
         q_act = gene_2;
         //From delta to gamma.
         aux += n2*zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         q_act = gene_1;
         //From gamma to delta.
         aux -= n2*zero_flow_matrix[genes[p_act]][genes[q_act]];
         change += aux*distAB;
      }
  }

  i_act = gene_1;
  j_act = gene_2;

  distAB = zero_dist_matrix[i_act][j_act];
  if(abs(distAB) > 0.000001){
      aux = 0;
      for( p = 0; p < n2; p++){
         p_act = it[p];
         q_act = gene_1;
         //From delta to gamma.
         aux += n2*zero_flow_matrix[genes[p_act]][genes[q_act]];
         q_act = gene_2;
         //From gamma to delta.
         aux -= n2*zero_flow_matrix[genes[p_act]][genes[q_act]];
      }
      p_act = j_act;
      q_act = gene_1;
      //From beta to alpha.
      aux += n24*zero_flow_matrix[genes[p_act]][genes[q_act]];
      for( q = 0; q < n2; q++){
         q_act = it[q];
         p_act = gene_1;
         //From gamma to delta.
         aux -= n2*zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         //From delta to gamma.
         aux += n2*zero_flow_matrix[genes[p_act]][genes[q_act]];
      }
      q_act = j_act;
      p_act = gene_1;
      //From alpha to beta.
      aux -= n24*zero_flow_matrix[genes[p_act]][genes[q_act]];
      change += aux*distAB;
   }

   i_act = gene_2;
   j_act = gene_1;
   distAB = zero_dist_matrix[i_act][j_act];
   if(abs(distAB) > 0.000001){
       aux = 0;
       for( p = 0; p < n2; p++){
          p_act = it[p];
          q_act = gene_1;
          //From gamma to delta.
          aux -= n2*zero_flow_matrix[genes[p_act]][genes[q_act]];
          q_act = gene_2;
          //From delta to gamma.
          aux += n2*zero_flow_matrix[genes[p_act]][genes[q_act]];
       }
       p_act = j_act;
       q_act = gene_2;
       //From beta to alpha.
       aux += n24*zero_flow_matrix[genes[p_act]][genes[q_act]];
       for( q = 0; q < n2; q++){
          q_act = it[q];
          p_act = gene_1;
          //From delta to gamma.
          aux += n2*zero_flow_matrix[genes[p_act]][genes[q_act]];
          p_act = gene_2;
          //From gamma to delta.
          aux -= n2*zero_flow_matrix[genes[p_act]][genes[q_act]];
       }
       q_act = j_act;
       p_act = gene_2;
       //From alpha to beta.
       aux -= n24*zero_flow_matrix[genes[p_act]][genes[q_act]];
       change += aux*distAB;
    }
    fitness += (long double)change/(n*n2);

    change = 0;
    aux = 0;
    i_act = gene_1;
    distAB = zero_dist_matrix[i_act][i_act];
    p_act = gene_1;
    aux -= zero_flow_matrix[genes[p_act]][genes[p_act]];
    p_act = gene_2;
    aux += zero_flow_matrix[genes[p_act]][genes[p_act]];
    change += aux*distAB;

    aux = 0;
    i_act = gene_2;
    distAB = zero_dist_matrix[i_act][i_act];
    p_act = gene_1;
    aux += zero_flow_matrix[genes[p_act]][genes[p_act]];
    p_act = gene_2;
    aux -= zero_flow_matrix[genes[p_act]][genes[p_act]];
    change += aux*distAB;

    fitness += (long double)change;

    free(it);
    return (fitness);
}

//Jointly computes f, f^1, f^2 and f^3 in order to save time. This function can only be used when the instance is asymmetric.
void evaluateQAP_functions123_asym(int *genes, long double *fitness){
      int i,j,p,q;
      int m_size = PermuSize;
      long int result1=0;
      long int result2=0;
      long int aux1_1,aux1_2,aux2_1,aux2_2;
      int n3=m_size-3;
      int n1=1-m_size;
      int n23=2*m_size-3;
      int n2=m_size-2;
      for (i=0;i<m_size;i++){
          for (j=0;j<i;j++){
              if(abs(zero_dist_matrix[i][j])<0.000001 && abs(zero_dist_matrix[j][i])<0.000001) continue;
              //min_plus=min+1;
              //max_plus=max+1;
              aux1_1=-m_pq[genes[i]][genes[j]];
              aux2_1=m_pq[genes[i]][genes[j]];

              aux1_2 = aux1_1;
              aux2_2 = aux2_1;

              for(q=0;q<genes[j];q++){
                  aux1_1 += -2*zero_flow_matrix[genes[i]][q];
              }
              for(q=genes[j]+1;q<m_size;q++){
                  aux1_1 += -2*zero_flow_matrix[genes[i]][q];
              }
              for(p=0;p<genes[i];p++){
                  aux1_1 += -2*zero_flow_matrix[p][genes[j]];
              }
              for(p=genes[i]+1;p<m_size;p++){
                  aux1_1 += -2*zero_flow_matrix[p][genes[j]];
              }

              for(q=0;q<genes[i];q++){
                  aux1_2 += -2*zero_flow_matrix[genes[j]][q];
              }
              for(q=genes[i]+1;q<m_size;q++){
                  aux1_2 += -2*zero_flow_matrix[genes[j]][q];
              }
              for(p=0;p<genes[j];p++){
                  aux1_2 += -2*zero_flow_matrix[p][genes[i]];
              }
              for(p=genes[j]+1;p<m_size;p++){
                  aux1_2 += -2*zero_flow_matrix[p][genes[i]];
              }

              //alpha case
              aux1_1+=n3*zero_flow_matrix[genes[i]][genes[j]];
              aux2_1+=n3*zero_flow_matrix[genes[i]][genes[j]];

              //beta case.
              aux1_1+=n1*zero_flow_matrix[genes[j]][genes[i]];
              aux2_1+=n3*zero_flow_matrix[genes[j]][genes[i]];

              //alpha case
              aux1_2+=n3*zero_flow_matrix[genes[j]][genes[i]];
              aux2_2+=n3*zero_flow_matrix[genes[j]][genes[i]];

              //beta case.
              aux1_2+=n1*zero_flow_matrix[genes[i]][genes[j]];
              aux2_2+=n3*zero_flow_matrix[genes[i]][genes[j]];

              //we multiply with distance here because the distance is a common factor in the calculation.
              result1+=aux1_1*zero_dist_matrix[i][j]+aux1_2*zero_dist_matrix[j][i];
              result2+=aux2_1*zero_dist_matrix[i][j]+aux2_2*zero_dist_matrix[j][i];
          }
      }


      fitness[0] = (long double)result1/(2*m_size);
      fitness[1] = (long double)result2/(2*(m_size-2));
      fitness[3] = (long double)evaluateQAP_function0(genes);
      fitness[2] = (long double)fitness[3]-(fitness[0]+fitness[1]);
}

//Jointly computes f, f^1, f^2 and f^3 in order to save time. This function can only be used when the instance is semi-symmetric (symmetric distance matrix).
void evaluateQAP_functions123_semi_sym1(int *genes, long double *fitness){
      int i,j,p,q;
      int m_size = PermuSize;
      long int result1=0;
      long int result2=0;
      long int aux1,aux2;
      int n3=m_size-3;
      int n1=1-m_size;
      int n23=2*m_size-3;
      int n2=m_size-2;
      for (i=0;i<m_size;i++){
          for (j=0;j<i;j++){
              if(abs(zero_dist_matrix[i][j])<0.000001) continue;
              //min_plus=min+1;
              //max_plus=max+1;
              aux1=-2*m_pq[genes[i]][genes[j]];
              aux2=2*m_pq[genes[i]][genes[j]];

              for(q=0;q<genes[j];q++){
                  aux1 += -2*zero_flow_matrix[genes[i]][q]-2*zero_flow_matrix[q][genes[i]];
              }
              for(q=genes[j]+1;q<m_size;q++){
                  aux1 += -2*zero_flow_matrix[genes[i]][q]-2*zero_flow_matrix[q][genes[i]];
              }
              for(p=0;p<genes[i];p++){
                  aux1 += -2*zero_flow_matrix[p][genes[j]]-2*zero_flow_matrix[genes[j]][p];
              }
              for(p=genes[i]+1;p<m_size;p++){
                  aux1 += -2*zero_flow_matrix[p][genes[j]]-2*zero_flow_matrix[genes[j]][p];
              }

              //alpha case
              aux1+=(n3+n1)*zero_flow_matrix[genes[i]][genes[j]];
              aux2+=2*n3*zero_flow_matrix[genes[i]][genes[j]];

              //beta case.
              aux1+=(n3+n1)*zero_flow_matrix[genes[j]][genes[i]];
              aux2+=2*n3*zero_flow_matrix[genes[j]][genes[i]];

              //we multiply with distance here because the distance is a common factor in the calculation.
              result1+=aux1*zero_dist_matrix[i][j];
              result2+=aux2*zero_dist_matrix[i][j];
          }
      }

      fitness[0] = (long double)result1/(2*m_size);
      fitness[1] = (long double)result2/(2*(m_size-2));
      fitness[3] = (long double)evaluateQAP_function0(genes);
      fitness[2] = (long double)fitness[3]-(fitness[0]+fitness[1]);
}

//Jointly computes f, f^1, f^2 and f^3 in order to save time. This function can only be used when the instance is semi-symmetric (symmetric flow matrix).
void evaluateQAP_functions123_semi_sym2(int *genes, long double *fitness){
      int i,j,p,q;
      int m_size = PermuSize;
      long int result1=0;
      long int result2=0;
      long int aux1_1,aux1_2,aux2_1,aux2_2;
      int n3=m_size-3;
      int n1=1-m_size;
      int n23=2*m_size-3;
      int n2=m_size-2;
      for (i=0;i<m_size;i++){
          for (j=0;j<i;j++){
              if(abs(zero_dist_matrix[i][j])<0.000001 && abs(zero_dist_matrix[j][i])<0.000001) continue;
              //min_plus=min+1;
              //max_plus=max+1;
              aux1_1=-m_pq[genes[i]][genes[j]];
              aux2_1=m_pq[genes[i]][genes[j]];

              aux1_2 = aux1_1;
              aux2_2 = aux2_1;

              for(q=0;q<genes[j];q++){
                  aux1_1 += -2*zero_flow_matrix[genes[i]][q];
              }
              for(q=genes[j]+1;q<m_size;q++){
                  aux1_1 += -2*zero_flow_matrix[genes[i]][q];
              }
              for(p=0;p<genes[i];p++){
                  aux1_1 += -2*zero_flow_matrix[p][genes[j]];
              }
              for(p=genes[i]+1;p<m_size;p++){
                  aux1_1 += -2*zero_flow_matrix[p][genes[j]];
              }

              for(q=0;q<genes[i];q++){
                  aux1_2 += -2*zero_flow_matrix[genes[j]][q];
              }
              for(q=genes[i]+1;q<m_size;q++){
                  aux1_2 += -2*zero_flow_matrix[genes[j]][q];
              }
              for(p=0;p<genes[j];p++){
                  aux1_2 += -2*zero_flow_matrix[p][genes[i]];
              }
              for(p=genes[j]+1;p<m_size;p++){
                  aux1_2 += -2*zero_flow_matrix[p][genes[i]];
              }

              //alpha case
              aux1_1+=(n3+n1)*zero_flow_matrix[genes[i]][genes[j]];
              aux2_1+=2*n3*zero_flow_matrix[genes[i]][genes[j]];

              //alpha case
              aux1_2+=(n3+n1)*zero_flow_matrix[genes[j]][genes[i]];
              aux2_2+=2*n3*zero_flow_matrix[genes[j]][genes[i]];

              //we multiply with distance here because the distance is a common factor in the calculation.
              result1+=aux1_1*zero_dist_matrix[i][j]+aux1_2*zero_dist_matrix[j][i];
              result2+=aux2_1*zero_dist_matrix[i][j]+aux2_2*zero_dist_matrix[j][i];
          }
      }

      fitness[0] = (long double)result1/(2*m_size);
      fitness[1] = (long double)result2/(2*(m_size-2));
      fitness[3] = (long double)evaluateQAP_function0(genes);
      fitness[2] = (long double)fitness[3]-(fitness[0]+fitness[1]);
}

//Jointly computes f, f^1, f^2 and f^3 in order to save time. This function can only be used when the instance is symmetric.
void evaluateQAP_functions123_sym(int *genes, long double *fitness){
      int i,j,p,q;
      int m_size = PermuSize;
      long int result1=0;
      long int result2=0;
      long int aux1,aux2;
      int n3=m_size-3;
      int n1=1-m_size;
      int n23=2*m_size-3;
      int n2=m_size-2;
      for (i=0;i<m_size;i++){
          for (j=0;j<i;j++){
              if(abs(zero_dist_matrix[i][j])<0.000001) continue;
              //min_plus=min+1;
              //max_plus=max+1;
              aux1=-m_pq[genes[i]][genes[j]];
              aux2=m_pq[genes[i]][genes[j]];

              for(q=0;q<genes[j];q++){
                  aux1 += -2*zero_flow_matrix[genes[i]][q];
              }
              for(q=genes[j]+1;q<m_size;q++){
                  aux1 += -2*zero_flow_matrix[genes[i]][q];
              }
              for(p=0;p<genes[i];p++){
                  aux1 += -2*zero_flow_matrix[p][genes[j]];
              }
              for(p=genes[i]+1;p<m_size;p++){
                  aux1 += -2*zero_flow_matrix[p][genes[j]];
              }

              //alpha case
              aux1+=(n3+n1)*zero_flow_matrix[genes[i]][genes[j]];
              aux2+=2*n3*zero_flow_matrix[genes[i]][genes[j]];

              //we multiply with distance here because the distance is a common factor in the calculation.
              result1+=2*aux1*zero_dist_matrix[i][j];
              result2+=2*aux2*zero_dist_matrix[i][j];
          }
      }

      fitness[0] = (long double)result1/(2*m_size);
      fitness[1] = (long double)result2/(2*(m_size-2));
      fitness[3] = (long double)evaluateQAP_function0(genes);
      fitness[2] = (long double)fitness[3]-(fitness[0]+fitness[1]);
}

//Jointly computes the f, f^1, f^2 and f^3 value difference between two neighbor solutions in order to save time.
//This function can only be used when the instance is asymmetric.
void evaluateQAP_change_function123_asym(long double * fitness, long double * fitness_new, int *genes, int *genes_swap, int gene_1, int gene_2){
  int n, i, j, p, q,n1,n2,n24,n3;
  int *it;
  int i_act, j_act, p_act, q_act, distAB;
  long int change1,change2,aux1,aux2;

  n = PermuSize;
  n1 = 1-n;
  n2 = n-2;
  n24 = 2*n-4;
  n3 = n-3;
  change1 = 0;
  change2 = 0;

  it=its[gene_1][gene_2];

  //Alpha, beta, gamma, delta, epsilon.
  for( i = 0; i < n2;i++){
     i_act = it[i];
     j_act = gene_1;
     distAB = zero_dist_matrix[i_act][j_act];
     if(abs(distAB) > 0.000001){
         aux1 = sum[genes[gene_1]][genes[gene_2]][genes[i_act]][0];
         aux2 = sum[genes[gene_1]][genes[gene_2]][genes[i_act]][0];

         p_act = i_act;
         q_act = gene_1;
         //From alpha to gamma.
         aux1 += n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         aux2-= n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         q_act = gene_2;
         //From gamma to alpha.
         aux1 -= n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         aux2+= n3*zero_flow_matrix[genes[p_act]][genes[q_act]];

         aux1 += -sum[genes[gene_1]][genes[gene_2]][genes[i_act]][1];
         aux2 += sum[genes[gene_1]][genes[gene_2]][genes[i_act]][1];

         q_act = i_act;
         p_act = gene_1;
         //From beta to delta.
         aux1 -= n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         aux2 -= n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         //From delta to beta.
         aux1 += n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         aux2 += n3*zero_flow_matrix[genes[p_act]][genes[q_act]];

         p_act = gene_1;
         q_act = gene_2;
         //From delta to gamma.
         aux1 -= 2*zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         q_act = gene_1;
         //From gamma to delta.
         aux1 += 2*zero_flow_matrix[genes[p_act]][genes[q_act]];
         change1 += aux1*distAB;
         change2 += aux2*distAB;
     }
     j_act = gene_2;
     distAB = zero_dist_matrix[i_act][j_act];
     if(abs(distAB) > 0.000001){
         aux1 = -sum[genes[gene_1]][genes[gene_2]][genes[i_act]][0];
         aux2 = -sum[genes[gene_1]][genes[gene_2]][genes[i_act]][0];

         p_act = i_act;
         q_act = gene_1;
         //From gamma to alpha.
         aux1 -= n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         aux2 += n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         q_act = gene_2;
         //From alpha to gamma.
         aux1 += n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         aux2 -= n3*zero_flow_matrix[genes[p_act]][genes[q_act]];

         aux1 += sum[genes[gene_1]][genes[gene_2]][genes[i_act]][1];
         aux2 += -sum[genes[gene_1]][genes[gene_2]][genes[i_act]][1];

         q_act = i_act;
         p_act = gene_1;
         //From delta to beta.
         aux1 += n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         aux2 += n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         //From beta to delta.
         aux1 -= n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         aux2 -= n3*zero_flow_matrix[genes[p_act]][genes[q_act]];

         p_act = gene_1;
         q_act = gene_2;
         //From gamma to delta.
         aux1 += 2*zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         q_act = gene_1;
         //From delta to gamma.
         aux1 -= 2*zero_flow_matrix[genes[p_act]][genes[q_act]];
         change1 += aux1*distAB;
         change2 += aux2*distAB;
      }
  }

  for( j = 0; j < n2;j++){
     j_act = it[j];
     i_act = gene_1;
     distAB = zero_dist_matrix[i_act][j_act];
     if(abs(distAB) > 0.000001){
         aux1 = -sum[genes[gene_1]][genes[gene_2]][genes[j_act]][0];
         aux2 = sum[genes[gene_1]][genes[gene_2]][genes[j_act]][0];

         p_act = j_act;
         q_act = gene_1;
         //From beta to delta.
         aux1 -= n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         aux2 -= n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         q_act = gene_2;
         //From delta to beta.
         aux1 += n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         aux2 += n3*zero_flow_matrix[genes[p_act]][genes[q_act]];

         aux1 += sum[genes[gene_1]][genes[gene_2]][genes[j_act]][1];
         aux2 += sum[genes[gene_1]][genes[gene_2]][genes[j_act]][1];

         q_act = j_act;
         p_act = gene_1;
         //From alpha to gamma.
         aux1 += n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         aux2 -= n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         //From gamma to alpha.
         aux1 -= n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         aux2 += n3*zero_flow_matrix[genes[p_act]][genes[q_act]];

         p_act = gene_1;
         q_act = gene_2;
         //From gamma to delta.
         aux1 += 2*zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         q_act = gene_1;
         //From delta to gamma.
         aux1 -= 2*zero_flow_matrix[genes[p_act]][genes[q_act]];
         change1 += aux1*distAB;
         change2 += aux2*distAB;
     }
     i_act = gene_2;
     distAB = zero_dist_matrix[i_act][j_act];
     if(abs(distAB) > 0.000001){
         aux1 = sum[genes[gene_1]][genes[gene_2]][genes[j_act]][0];
         aux2 = -sum[genes[gene_1]][genes[gene_2]][genes[j_act]][0];

         p_act = j_act;
         q_act = gene_1;
         //From delta to beta.
         aux1 += n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         aux2 += n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         q_act = gene_2;
         //From beta to delta.
         aux1 -= n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         aux2 -= n3*zero_flow_matrix[genes[p_act]][genes[q_act]];

         aux1 += -sum[genes[gene_1]][genes[gene_2]][genes[j_act]][1];
         aux2 += -sum[genes[gene_1]][genes[gene_2]][genes[j_act]][1];

         q_act = j_act;
         p_act = gene_1;
         //From gamma to alpha.
         aux1 -= n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         aux2 += n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         //From alpha to gamma.
         aux1 += n1*zero_flow_matrix[genes[p_act]][genes[q_act]];
         aux2 -= n3*zero_flow_matrix[genes[p_act]][genes[q_act]];

         p_act = gene_1;
         q_act = gene_2;
         //From delta to gamma.
         aux1 -= 2*zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         q_act = gene_1;
         //From gamma to delta.
         aux1 += 2*zero_flow_matrix[genes[p_act]][genes[q_act]];
         change1 += aux1*distAB;
         change2 += aux2*distAB;
      }
  }

  i_act = gene_1;
  j_act = gene_2;

  distAB = zero_dist_matrix[i_act][j_act];
  if(abs(distAB) > 0.000001){
      aux1 = -2*sum_exp[genes[gene_1]][genes[gene_2]][0];

      p_act = j_act;
      q_act = gene_1;
      //From beta to alpha.
      aux1 += n24*zero_flow_matrix[genes[p_act]][genes[q_act]];

      aux1 += 2*sum_exp[genes[gene_1]][genes[gene_2]][1];

      q_act = j_act;
      p_act = gene_1;
      //From alpha to beta.
      aux1 -= n24*zero_flow_matrix[genes[p_act]][genes[q_act]];
      change1 += aux1*distAB;
   }

   i_act = gene_2;
   j_act = gene_1;
   distAB = zero_dist_matrix[i_act][j_act];
   if(abs(distAB) > 0.000001){
       aux1 = 2*sum_exp[genes[gene_1]][genes[gene_2]][0];

       p_act = j_act;
       q_act = gene_2;
       //From beta to alpha.
       aux1 += n24*zero_flow_matrix[genes[p_act]][genes[q_act]];

       aux1 += -2*sum_exp[genes[gene_1]][genes[gene_2]][1];

       q_act = j_act;
       p_act = gene_2;
       //From alpha to beta.
       aux1 -= n24*zero_flow_matrix[genes[p_act]][genes[q_act]];
       change1 += aux1*distAB;
    }

    fitness_new[0] = fitness[0] + (long double)change1/(2*n);
    fitness_new[1] = fitness[1] + (long double)change2/(2*n2);
    fitness_new[3] = (long double)evaluateQAP_change_function0(fitness[3],gene_1,gene_2);
    fitness_new[2] = (long double)fitness_new[3]-(fitness_new[0]+fitness_new[1]);
}

//Jointly computes the f, f^1, f^2 and f^3 value difference between two neighbor solutions in order to save time.
//This function can only be used when the instance is semi-symmetric (symmetric distance matrix).
void evaluateQAP_change_function123_semi_sym1(long double * fitness, long double * fitness_new, int *genes, int *genes_swap, int gene_1, int gene_2){
  int n, i, j, p, q,n1,n2,n24,n3;
  int *it;
  int i_act, j_act, p_act, q_act, distAB;
  long int change2,aux2;

  n = PermuSize;
  n1 = 1-n;
  n2 = n-2;
  n24 = 2*n-4;
  n3 = n-3;
  change2 = 0;

  it=its[gene_1][gene_2];

  //Alpha, beta, gamma, delta, epsilon.
  for( i = 0; i < n2;i++){
     i_act = it[i];
     j_act = gene_1;
     distAB = zero_dist_matrix[i_act][j_act];
     if(abs(distAB) > 0.000001){
         aux2 = 2*sum[genes[gene_1]][genes[gene_2]][genes[i_act]][0];

         p_act = i_act;
         q_act = gene_1;
         //From alpha to gamma.
         aux2-= 2*n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         q_act = gene_2;
         //From gamma to alpha.
         aux2+= 2*n3*zero_flow_matrix[genes[p_act]][genes[q_act]];

         aux2 += 2*sum[genes[gene_1]][genes[gene_2]][genes[i_act]][1];

         q_act = i_act;
         p_act = gene_1;
         //From beta to delta.
         aux2 -= 2*n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         //From delta to beta.
         aux2 += 2*n3*zero_flow_matrix[genes[p_act]][genes[q_act]];

         change2 += aux2*distAB;
     }
     j_act = gene_2;
     distAB = zero_dist_matrix[i_act][j_act];
     if(abs(distAB) > 0.000001){ //CONTINUAR (Dist simetrico)
         aux2 = -2*sum[genes[gene_1]][genes[gene_2]][genes[i_act]][0];

         p_act = i_act;
         q_act = gene_1;
         //From gamma to alpha.
         aux2 += 2*n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         q_act = gene_2;
         //From alpha to gamma.
         aux2 -= 2*n3*zero_flow_matrix[genes[p_act]][genes[q_act]];

         aux2 += -2*sum[genes[gene_1]][genes[gene_2]][genes[i_act]][1];

         q_act = i_act;
         p_act = gene_1;
         //From delta to beta.
         aux2 += 2*n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         p_act = gene_2;
         //From beta to delta.
         aux2 -= 2*n3*zero_flow_matrix[genes[p_act]][genes[q_act]];

         change2 += aux2*distAB;
      }
  }

    fitness_new[0] = fitness[0];
    fitness_new[1] = fitness[1] + (long double)change2/(2*n2);
    fitness_new[3] = (long double)evaluateQAP_change_function0(fitness[3],gene_1,gene_2);
    fitness_new[2] = (long double)fitness_new[3]-(fitness_new[0]+fitness_new[1]);

}

//Jointly computes the f, f^1, f^2 and f^3 value difference between two neighbor solutions in order to save time.
//This function can only be used when the instance is semi-symmetric (symmetric flow matrix).
void evaluateQAP_change_function123_semi_sym2(long double * fitness, long double * fitness_new, int *genes, int *genes_swap, int gene_1, int gene_2){
  int n, i, j, p, q,n1,n2,n24,n3;
  int *it;
  int i_act, j_act, p_act, q_act, distAB;
  long int change2,aux2;

  n = PermuSize;
  n1 = 1-n;
  n2 = n-2;
  n24 = 2*n-4;
  n3 = n-3;
  change2 = 0;

  it=its[gene_1][gene_2];

  //Alpha, beta, gamma, delta, epsilon.
  for( i = 0; i < n2;i++){
     i_act = it[i];
     j_act = gene_1;
     distAB = zero_dist_matrix[i_act][j_act];
     if(abs(distAB) > 0.000001){
         aux2 = 2*sum[genes[gene_1]][genes[gene_2]][genes[i_act]][0];

         p_act = i_act;
         q_act = gene_1;
         //From alpha to gamma.
         aux2-= 2*n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         q_act = gene_2;
         //From gamma to alpha.
         aux2+= 2*n3*zero_flow_matrix[genes[p_act]][genes[q_act]];

         change2 += aux2*distAB;
     }
     j_act = gene_2;
     distAB = zero_dist_matrix[i_act][j_act];
     if(abs(distAB) > 0.000001){
         aux2 = -2*sum[genes[gene_1]][genes[gene_2]][genes[i_act]][0];

         p_act = i_act;
         q_act = gene_1;
         //From gamma to alpha.
         aux2 += 2*n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         q_act = gene_2;
         //From alpha to gamma.
         aux2 -= 2*n3*zero_flow_matrix[genes[p_act]][genes[q_act]];

         change2 += aux2*distAB;
      }
  }

  for( j = 0; j < n2;j++){
     j_act = it[j];
     i_act = gene_1;
     distAB = zero_dist_matrix[i_act][j_act];
     if(abs(distAB) > 0.000001){
         aux2 = 2*sum[genes[gene_1]][genes[gene_2]][genes[j_act]][0];

         p_act = j_act;
         q_act = gene_1;
         //From beta to delta.
         aux2 -= 2*n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         q_act = gene_2;
         //From delta to beta.
         aux2 += 2*n3*zero_flow_matrix[genes[p_act]][genes[q_act]];

         change2 += aux2*distAB;
     }
     i_act = gene_2;
     distAB = zero_dist_matrix[i_act][j_act];
     if(abs(distAB) > 0.000001){
         aux2 = -2*sum[genes[gene_1]][genes[gene_2]][genes[j_act]][0];

         p_act = j_act;
         q_act = gene_1;
         //From delta to beta.
         aux2 += 2*n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         q_act = gene_2;
         //From beta to delta.
         aux2 -= 2*n3*zero_flow_matrix[genes[p_act]][genes[q_act]];

         change2 += aux2*distAB;
      }
  }

    fitness_new[0] = fitness[0];
    fitness_new[1] = fitness[1] + (long double)change2/(2*n2);
    fitness_new[3] = (long double)evaluateQAP_change_function0(fitness[3],gene_1,gene_2);
    fitness_new[2] = (long double)fitness_new[3]-(fitness_new[0]+fitness_new[1]);
}

//Jointly computes the f, f^1, f^2 and f^3 value difference between two neighbor solutions in order to save time.
//This function can only be used when the instance is symmetric.
void evaluateQAP_change_function123_sym(long double * fitness, long double * fitness_new, int *genes, int *genes_swap, int gene_1, int gene_2){
  int n, i, j, p, q,n1,n2,n24,n3;
  int *it;
  int i_act, j_act, p_act, q_act, distAB;
  long int change2,aux2;

  n = PermuSize;
  n1 = 1-n;
  n2 = n-2;
  n24 = 2*n-4;
  n3 = n-3;
  change2 = 0;

  it=its[gene_1][gene_2];

  //Alpha, beta, gamma, delta, epsilon.
  for( i = 0; i < n2;i++){
     i_act = it[i];
     j_act = gene_1;
     distAB = zero_dist_matrix[i_act][j_act];
     if(abs(distAB) > 0.000001){
         aux2 = 4*sum[genes[gene_1]][genes[gene_2]][genes[i_act]][0];

         p_act = i_act;
         q_act = gene_1;
         //From alpha to gamma.
         aux2 -= 4*n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         q_act = gene_2;
         //From gamma to alpha.
         aux2 += 4*n3*zero_flow_matrix[genes[p_act]][genes[q_act]];

         change2 += aux2*distAB;
     }
     j_act = gene_2;
     distAB = zero_dist_matrix[i_act][j_act];
     if(abs(distAB) > 0.000001){
         aux2 = -4*sum[genes[gene_1]][genes[gene_2]][genes[i_act]][0];

         p_act = i_act;
         q_act = gene_1;
         //From gamma to alpha.
         aux2 += 4*n3*zero_flow_matrix[genes[p_act]][genes[q_act]];
         q_act = gene_2;
         //From alpha to gamma.
         aux2 -= 4*n3*zero_flow_matrix[genes[p_act]][genes[q_act]];

         change2 += aux2*distAB;
      }
  }

    fitness_new[0] = fitness[0];
    fitness_new[1] = fitness[1] + (long double)change2/(2*n2);
    fitness_new[3] = (long double)evaluateQAP_change_function0(fitness[3],gene_1,gene_2);
    fitness_new[2] = (long double)fitness_new[3]-(fitness_new[0]+fitness_new[1]);
}

//Auxiliary function that can be called outside the current file. It efficiently computes f, f^1, f^2 or f^3.
long double EvaluateQAP(int fitnessFunctionId, int * genes){
    long double fitness = INT_MAX;
    switch (fitnessFunctionId){
        case 1:
            fitness = evaluateQAP_function1_opt(genes);
            break;
        case 2:
            fitness = evaluateQAP_function2_opt(genes);
            break;
        case 3:
            fitness = evaluateQAP_function3_opt(genes);
            break;
        default:
            fitness = evaluateQAP_function0(genes);
            break;
    }

    return fitness;
}

//Auxiliary function that can be called outside the current file. It jointly computes f, f^1, f^2 and f^3.
void EvaluateQAP_all(long double * fitness, int * genes){
    if(Symmetry == 0){
        evaluateQAP_functions123_asym(genes, fitness);
    }else if(Symmetry == 1){
        evaluateQAP_functions123_semi_sym1(genes, fitness);
    }else if(Symmetry == 2){
        evaluateQAP_functions123_semi_sym2(genes, fitness);
    }else{
        evaluateQAP_functions123_sym(genes, fitness);
    }
}

//Auxiliary function that can be called outside the current file. It computes the f, f^1, f^2 or f^3 value difference between two neighbor solutions.
long double EvaluateQAP_change(int fitnessFunctionId, long double fitness, int * genes, int * genes_swap, int gene_1, int gene_2){
    switch (fitnessFunctionId){
        case 1:
            fitness = evaluateQAP_change_function1_opt(fitness, genes, genes_swap, gene_1, gene_2);
            break;
        case 2:
            fitness = evaluateQAP_change_function2_opt(fitness, genes, genes_swap, gene_1, gene_2);
            break;
        case 3:
            fitness = evaluateQAP_change_function3_opt(fitness, genes, genes_swap, gene_1, gene_2);
            break;
        default:
            fitness = evaluateQAP_change_function0(fitness, gene_1, gene_2);
            break;
    }

    return fitness;
}

//Auxiliary function that can be called outside the current file. It jointly computes the f, f^1, f^2 and f^3 value difference between two neighbor solutions.
void EvaluateQAP_change_all(long double * fitness, long double * fitness_new, int * genes, int * genes_swap, int gene_1, int gene_2){
    if(Symmetry == 0){
        evaluateQAP_change_function123_asym(fitness, fitness_new, genes, genes_swap, gene_1, gene_2);
    }else if(Symmetry == 1){
        evaluateQAP_change_function123_semi_sym1(fitness, fitness_new, genes, genes_swap, gene_1, gene_2);
    }else if(Symmetry == 2){
        evaluateQAP_change_function123_semi_sym2(fitness, fitness_new, genes, genes_swap, gene_1, gene_2);
    }else{
        evaluateQAP_change_function123_sym(fitness, fitness_new, genes, genes_swap, gene_1, gene_2);
    }
}
