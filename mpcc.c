#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#define MAX_FLOWS 10
#define MAX_SUBFLOWS 10
#define K 10
#define L 1

typedef struct subflow;

typedef struct {
  double capacity;
  struct subflow** subs;
  int subs_count;
  double loss;
} link;

typedef struct {
  link** path;
  int path_len;
  double rate;
  double loss;
} subflow;

typedef struct {
  subflow* subs;
  int subs_count;

  double * x[K];
  double * v[K];
  double * t;
  int probing;
  int phase;
}flow;

void compute_link_loss(link* l){
  int i,j;
  double d = 0;
  
  for (j=0;j<l->subs_count;j++)
    d += ((subflow*)l->subs[j])->rate;
  
  if (d==0||d<l->capacity)
    l->loss = 0;
  else 
    l->loss = ((double)d-l->capacity)/d;

  //  printf("[L%p %f]",l,l->loss);
};

void compute_subflow_loss(subflow* s){
  int i,j;
  double d;

  d = 0;
  //iterate all links of subflow s
  for (j = 0;j<s->path_len;j++){
    d += s->path[j]->loss;
  }

  s->loss = d;
}

double utility(flow* f){
  double u1 = 0, u2 = 0;
  int i;
  
  //  printf("U(");

  for (i = 0;i<f->subs_count;i++){
    //printf("%f,%f\t",f->subs[i].rate,f->subs[i].loss);
    u1 += f->subs[i].rate;
    u2 += f->subs[i].rate * (pow(1+f->subs[i].loss,5)-1);
  }
  //  printf(")=%f",u1-u2);

  return u1-u2;
}

void run_cc(flow* f){
  int i,j;
  double delta = 10,epsilon = 0.05;

  if(f->probing){
    if (f->phase>0){
      double u = utility(f);

      //printf("\nGRAD (");
      for (i = 0;i<f->subs_count;i++){
	//compute gradient for previous phase
	f->x[f->phase-1][i] = (double)f->subs_count / delta * u * f->v[f->phase-1][i];
	
	//printf("%f ",f->x[f->phase-1][i]);
      
      //correct rate for last phase
	f->subs[i].rate = f->t[i];
      }
      //printf(")\n");

      if (f->phase==K){
	//compute average gradient
	for (i=1;i < K;i++){
	  for (j=0;j<f->subs_count;j++)
	    f->x[0][j] += f->x[i][j];
	}

	//printf("\nCHOOSE GRAD_AVG(\t");
	for (j=0;j<f->subs_count;j++){
	  f->x[0][j] /= K;
	  //printf("%f ",f->x[0][j]);

	  //adapt the rate
	  if (f->subs[j].rate + epsilon * f->x[0][j]>1)
	    f->subs[j].rate += epsilon * f->x[0][j];
	  else
	    f->subs[j].rate = 1;
	}

	f->probing = 0;
	f->phase = 0;//rand()%L;
	return;
      }
    }
    
    //choose random unit vector
    double norm = 0;

    for (i = 0;i<f->subs_count;i++){
      f->v[f->phase][i] = (double)rand()-RAND_MAX/2;
      norm += f->v[f->phase][i]*f->v[f->phase][i];
    }
    norm = sqrt(norm);

    //printf("\nTRY (");
    for (i = 0;i<f->subs_count;i++){
      f->v[f->phase][i] /= norm;
      //printf("%f ",f->v[f->phase][i]);
    }
    //printf(")");

    //now update rate to reflect this vector, try it in the next cuanta

    for (i = 0;i<f->subs_count;i++){
      if (f->subs[i].rate > -delta * f->v[f->phase][i]){
	f->subs[i].rate = f->subs[i].rate + delta * f->v[f->phase][i];    
      }
      else f->subs[i].rate = 1;
    }


    //increase phase
    f->phase++;
  }
  else {//not probing
    f->phase++;
    if (f->phase>L){
      for (i=0;i<f->subs_count;i++)
	f->t[i] = f->subs[i].rate;
      f->probing = 1;
      f->phase = 0;//rand()%K;
    }
  }
}

void add_link_to_sub(subflow* s, link* l){
  s->path[s->path_len++] = l;
  l->subs[l->subs_count++] = (struct subflow*)s;
}

int main(int argc,char** argv){
  double delta = 1,epsilon =1;
  int i,j;
  
  //create topology;

  link* links[3];
  flow* flows[3];

  for (i=0;i<3;i++){
    links[i] = (link*) calloc(1,sizeof(link));
    links[i]->capacity = 100;
    links[i]->subs = (struct subflow**)calloc(3,sizeof(subflow*));
  }

  for (i=0;i<3;i++){
    flows[i] = (flow*) calloc(1,sizeof(flow));

    flows[i]->subs_count = 2;

    for (j=0;j<K;j++){
      flows[i]->x[j] = (double*)calloc(flows[i]->subs_count,sizeof(double));
      flows[i]->v[j] = (double*)calloc(flows[i]->subs_count,sizeof(double));
    }
    flows[i]->t = (double*)calloc(flows[i]->subs_count,sizeof(double));
    
    flows[i]->subs = (subflow*)calloc(flows[i]->subs_count,sizeof(subflow));
    
    flows[i]->subs[0].path = (link**)calloc(1,sizeof(link*));
    add_link_to_sub(&(flows[i]->subs[0]),links[i]);
    flows[i]->subs[0].rate = 110;

    if (flows[i]->subs_count>1){
      flows[i]->subs[1].path = (link**)calloc(2,sizeof(link*));
      add_link_to_sub(&(flows[i]->subs[1]),links[(i+1)%3]);
      add_link_to_sub(&(flows[i]->subs[1]),links[(i+2)%3]);

      flows[i]->subs[1].rate = 110;
    }
  }

  //now start the simulator

  j = 100000;
  int k;
  while (j-->0){
    //compute link loss rates
    for (i=0;i<3;i++)
      //link loss 
      compute_link_loss(links[i]);
    
    for (i=0;i<3;i++){
      //compute subflow loss rates
      for (k=0;k<flows[i]->subs_count;k++)
	compute_subflow_loss(&(flows[i]->subs[k]));

      //apply congestion control
      run_cc(flows[i]);
    }

    //    if (j%10==0){
      for (i=0;i<3;i++)
	printf("[ %f %f ]\t",flows[i]->subs[0].rate,flows[i]->subs[1].rate);

      printf("\n");
      //};
  }
}
