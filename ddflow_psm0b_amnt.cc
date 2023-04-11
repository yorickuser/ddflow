// file ddflow_psmb0b.cc
// C++ program for simulation of DD-flows with PSM
// (Polymorphic stochastic model, Dieckmann and Law (1996))
// for sexual populations (see Ito and Dieckmann (2007)).
//
// Written by Hiroshi C. Ito (2021)
// Email: hiroshibeetle@gmail.com

#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

/*_/_/_/_/_/_/_/_/_/ Model prameters _/_/_/_/_/_/_/_/*/
#define SEXUAL_REPRODUCTION 1
//1:sexual reproduction, 0:asexual reproduction

double gradw=4;  /* fitness gradient for fundamental trait (beta=gradw*max_w) */
double assortative_rate = 2.0;   /* curvature of mating probability function */
double self_reproduction_prob=0; /* probalility for self reproduction when finding mating partners is difficult*/

const int size_p=120; /* number of loci for niche trait */
const int size_w = 50;/* number of loci for fundamental trait */
const int size_m = 40;/* number of loci for mating trait (former half: male display trait, latter half: female preference trait)*/
const double max_w = 0.5;
const int amp_size_w = 4;

double delta_p = 1.0/double(size_p*2.0);      /* mutation step size for p */
double delta_w = max_w/double(size_w*2.0);  /* mutation step size for w */

/*cupling of loci on the same chromosome
value 1 corresponds to that all loci exists on different chromosomes */
int chro_p = 10;     
int chro_w = 10;
int chro_m = 10;

/* mutation probability at one locus */
double mutation_rate_p = 0.000025;
double mutation_rate_w = 0.00001;
double mutation_rate_m = 0.005;   

/* genotype of individuals at the initial condition */
/* number of allele 1 in the gene */
int first_p = size_p/2; 
int first_w = 0;

/* carrying capacity */
double R0=500; // peak
double mR=0.5; // peak position
double sdR=0.15; // sd

/* competition kernel */
double alpha=0.04; // sd

int seed = 13;  /* random seed */

/* species identification */
const int spe_number = 100;
int resol=2;
double wave_edge = 2;


int phenotype_interval = 1;
int count_interval=20000;

int Flag_output = 0;
int Flag_out_species =0;
int Flag_out_phenotype =0;
int Flag_stop = 1;

double offset_w=0.0;

const int g_size_x = size_p*2+1;
const int g_size_y = amp_size_w*(size_w*2+1);


/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/
/*_/_/ Functions for interaction and reproduction for individuals _/_/*/
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/


#include <cmath>
#include <iostream>
using namespace std;
const double Rand =  RAND_MAX;
#define PI 3.14159265358979323846


//_/_/_/_/_/_/_/_/_/_/_/_/ Definition of individual _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_

int random_int(int);

/*  Definition of a class for one individual (diploid) */
class Ind {
  public:
    Ind(int p[2*size_p],  int w[2*size_w],int m[2*size_m]);
    void reset(int p[2*size_p],int w[2*size_w],int m[2*size_m]);
    
    ~Ind();
    unsigned int alive :1;     /* 1:alive 0:dead */
    unsigned int sex;
    int gene_p[2*size_p];      /* Loci for niche trait. The former half corresponds to one chromosome, the latter half corresponds to the other chromosome. */ 
    int gene_w[2*size_w];     /* Loci for fundamental trait */
    int gene_m[2*size_m];     /* Loci for mating trait */
    double value_p;           /* niche trait */
    double value_w;           /* fundamental trait */
    double value_m;           /* mating trait */
    int v_p;                  /* ID for niche trait */
    double birth;             /* birth rate */
    double death;             /* death rate */
    
  } ;      
    
/* Constructer */
Ind::Ind(int p[2*size_p],int w[2*size_w],int m[2*size_m]) {

  alive = 1;
  sex=random_int(1);
  value_p = 0;
  value_w = delta_w;
  value_m = 0;
  birth=0;
  death=0;  
 
  for(int i=0; i<2*size_p; i++){
    gene_p[i] = p[i];
    value_p += delta_p*double(p[i]);
    v_p+=p[i];
  }
   for(int i=0; i<2*size_w; i++){gene_w[i] = w[i];value_w += delta_w*double(w[i]);}
 
  for(int i=0; i<2*size_m; i++)gene_m[i] = m[i];

    
  
}

/*  Resetting function for recycle of instances of dead individuals*/
void Ind::reset(int p[2*size_p],int w[2*size_w], int m[2*size_m]) {

  alive = 1;
  sex=random_int(1);
  value_p = 0;
  v_p=0;
  value_w = delta_w;
  value_m = 0;
  birth=0;
  death=0;

  for(int i=0; i<2*size_p; i++){
    gene_p[i] = p[i];
    value_p += delta_p*double(p[i]);
    v_p+=p[i];
  }
  for(int i=0; i<2*size_w; i++){gene_w[i] = w[i];value_w += delta_w*double(w[i]);}
  for(int i=0; i<2*size_m; i++)gene_m[i] = m[i];
 
}




int step=0;
double timen=0;
const int Z_division = size_p*2+1;  /* division of resource space for graphics */
const double delta_z = 1.0/double(Z_division); 

double Phe[size_p*2+1],Com[size_p*2+1],Com_sum;
double KK[size_p*2+1],comm[size_p*2+1];
int phenotype[size_p*2+1][amp_size_w*(size_w*2+1)];  /* number of individual in each phenotype */
double R[Z_division];                   /* Resource distribution */ 
double U[Z_division];                   /* Total utilization distribution */

const int MAX_NUMBER = 50000;             /* limitation of number of individual */
Ind *l[MAX_NUMBER];                  /* ID for individuals */

int number_end;                      /* total number of being used ID < MAX_NUMBER */
int number_end_buf;
int number;                          /* number of individuals */
int parent[2][MAX_NUMBER],parents[MAX_NUMBER];               /* ID for potential parents */
int baby_id=0;
int parent_buf[2][MAX_NUMBER];               /* ID for potential parents */

double vp, vw, vpw;
double gp,gl,gw,gp0,gp1,gw0,gw1,gm0x,gm1x,gm0y,gm1y,dist_m, gl_0, gl_1,maxw;

void number_count(void);
int prob(double); 
int prob2(double m_rate);
void Birth(Ind *ind0, Ind *ind1);
void random_sort(int size, int a[]);

void Clone(Ind *ind);
double f_Norm(double x, double m, double sd);

void birth(int gene_p[2*size_p],int gene_w[2*size_w],int  gene_m[2*size_m]);

int gene_p[2*size_p],gene_w[2*size_w],gene_m[2*size_m];

void mutation_only_stepwise(int genesize, double mutation_rate,int gene_parent[],int gene_baby[]);
void recombination_mutation_stepwise(int genesize, double mutation_rate,int chro_same,int gene_parent0[],int gene_parent1[],int gene_baby[]);

int mutation_stepwise(int x, double mutation_rate);

void Event(void);
int nth = 1;


//_/_/_/_/_/_/_/_/_/ Some primitive functions _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

/* random sorting of the array a[] */
void random_sort(int size, int a[]){
  int j,buf;
  int random_int(int);
  for(int i=0; i<size-1; i++){
    j = random_int(size-1-i);
    buf = a[j];
    a[j] = a[size-1-i];
    a[size-1-i] = buf;
  } 
}

/*  return a random integer ranging from 0 to M */
 int random_int(int M){
   return rand()%(M+1);
 }

/* return 1 at probability p, return -1 at probability (1-p) */
int prob(double p){
  if (rand()/Rand < p) return 1;
  else return -1;
}

/* return 0 at probability p, return 1 at probability 0.5*(1-p), return -1 at probability 0.5*(1-p) */
int prob2(double m_rate){
    double a=rand()/Rand;
    if (a > m_rate) return 0;
    else{
        if(a<0.5*m_rate)return 1;
        else return -1;
    }
}


/* Gausian function */
double f_Norm(double x, double m, double sd){
  return (1/(sd*sqrt(2*PI)))*exp(-(x-m)*(x-m)/(2*sd*sd));
}



//_/_/_/_/_/_/_/_/_/_/_/_/  Resouce Competition  _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
/* competition kernel */
inline double com(double x0, double x1){
    return exp(-1.0*(x0-x1)*(x0-x1)/(2*alpha*alpha));
}

/* resource distribution */
inline double K(double s0){
    return R0*exp(-1.0*(s0-mR)*(s0-mR)/(2*sdR*sdR));
}


//_/_/_/_/_/_/_/_/_/_/_/_/  Reproduction   _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

inline double mate_distance(Ind *ind0, Ind *ind1){
double dist;
  double value;
  double m0,m1;
  double y;
  
/* calculation of distance in mating traits and niche trait between two individuals */
  m0 = m1 =0;
  dist=0;
   
   for(int i=0; i<size_m/2; i++){
     m0 = (ind0->gene_m[i]+ind0->gene_m[i+size_m]);
     m1 = (ind1->gene_m[i+size_m/2]+ind1->gene_m[i+size_m/2+size_m]);
      dist += double((m0-m1)*(m0-m1));
  }
 
   //value = com(ind0->value_p,ind1->value_p)*exp(-1*dist/(assortative_rate*assortative_rate));
  value = comm[abs(ind0->v_p-ind1->v_p)]*exp(-1*dist/(2*assortative_rate*assortative_rate));
   return value;

}

int n_par[2];
int n_pars;
double sum_n,sum_b,sum_bm,sum_d;

inline void calc_npar(void){
  int sex;
  sum_n=sum_b=sum_d=0;
  n_pars=n_par[0]=n_par[1]= 0;

  for(int i=0; i<= number_end; i++){
    if(l[i]->alive > 0){
      l[i]->birth =1+gradw*(l[i]->value_w+offset_w);
      if(l[i]->birth<0)l[i]->birth=0;
      sex=l[i]->sex;
            
      parents[n_pars]=i;
      n_pars++;
      parent[sex][n_par[sex]] = i; 
      n_par[sex]++;
      if(sex==0)sum_n+=l[i]->death;
      else sum_n+=2.0*l[i]->birth+l[i]->death;
      
      sum_d+=l[i]->death;
      if(l[i]->sex==1)sum_b+=l[i]->birth;
      
    }
  }
}

inline void do_death(void){
  int k;
  double buf,rn;
  k=-1;
  rn=rand()/Rand;
  buf=0;
  while(rn>buf){
    k++;            
    buf+=l[parents[k]]->death/sum_d;
  }
  
  k=parents[k];
  
  l[k]->alive=0;
  for(int i=0;i<n_pars;i++){
    //l[parents[i]]->death-=com(l[k]->value_p,l[parents[i]]->value_p)/K(l[parents[i]]->value_p);
    l[parents[i]]->death-=comm[abs(l[k]->v_p-l[parents[i]]->v_p)]/KK[l[parents[i]]->v_p];
  }
  
}

inline void do_birth_sexual(int pf){
  int pm;
  double buf,rn;
  
  sum_bm=0;
  for(int i=0;i<n_par[0];i++){
    //double buff=(l[parent[0][i]]->birth)*mate_distance(l[parent[0][i]],l[pf]);
    l[parent[0][i]]->birth=(l[parent[0][i]]->birth)*mate_distance(l[parent[0][i]],l[pf]);
    sum_bm+=l[parent[0][i]]->birth;
  }
  if(prob(self_reproduction_prob/(sum_bm+self_reproduction_prob))==1){
    Clone(l[pf]);
  }
  else{
    
    pm=-1;
    rn=rand()/Rand;
    buf=0;    
    while(rn>buf){
      pm++;
      buf+=l[parent[0][pm]]->birth/sum_bm;      
    }
    pm=parent[0][pm];
    Birth(l[pm],l[pf]);
  }

}

void  Event(void) {

    int n_par_buf[2], m_count,pf,pm,k,kn,n_par0[2];
    int id,done;
    double bsize,m_ind,rn,buf;
    number_end_buf = number_end; 
    
    calc_npar();

    timen+=-(1.0/sum_n)*log((rand()+1.0)/(Rand+2.0));
    
    rn=rand()/Rand;
    if(rn<sum_d/sum_n){
      do_death();	
    }
    else{ 
        pf=-1;
        rn=rand()/Rand;
        buf=0;
        while(rn>buf){
            pf++;            
            buf+=(l[parent[1][pf]]->birth)/sum_b;
            
        }
        pf=parent[1][pf];
        
        if(SEXUAL_REPRODUCTION){

	  do_birth_sexual(pf);
        }
        else Clone(l[pf]);
        

        for(int i=0;i<n_pars;i++){
	  //l[parents[i]]->death+=com(l[baby_id]->value_p,l[parents[i]]->value_p)/K(l[parents[i]]->value_p);
	  l[parents[i]]->death+=comm[abs(l[baby_id]->v_p-l[parents[i]]->v_p)]/KK[l[parents[i]]->v_p];
        }
	
        for(int i=0;i<n_pars;i++){
	  //l[baby_id]->death+=com(l[baby_id]->value_p,l[parents[i]]->value_p);
	  l[baby_id]->death+=comm[abs(l[baby_id]->v_p-l[parents[i]]->v_p)];
        }
        l[baby_id]->death/=KK[l[baby_id]->v_p];


    number_end = number_end_buf; 
    }


    
}




//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

/* function for clonal reproduction */
void Clone (Ind *ind){
    mutation_only_stepwise(size_p, mutation_rate_p, ind->gene_p, gene_p);
    mutation_only_stepwise(size_w, mutation_rate_w, ind->gene_w, gene_w);
    mutation_only_stepwise(size_m, mutation_rate_m, ind->gene_m, gene_m);
    birth(gene_p,gene_w,gene_m);
}

inline void mutation_only_stepwise(int genesize, double mutation_rate,int gene_parent[],int gene_baby[]){
  for(int i=0; i<2*genesize; i++){
    gene_baby[i] = mutation_stepwise(gene_parent[i],mutation_rate);
  }
}

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

/*  function for sexual reproduction  */
void Birth(Ind *ind0, Ind *ind1){
    int c0, c1;
    recombination_mutation_stepwise(size_p,mutation_rate_p,chro_p,ind0->gene_p ,ind1->gene_p, gene_p);
    recombination_mutation_stepwise(size_w,mutation_rate_w,chro_w,ind0->gene_w ,ind1->gene_w, gene_w);
    recombination_mutation_stepwise(size_m,mutation_rate_m,chro_m,ind0->gene_m ,ind1->gene_m, gene_m);

    birth(gene_p,gene_w,gene_m);
}


/* function for recombination and mutation (stepwise mutation)*/
inline void recombination_mutation_stepwise(int genesize, double mutation_rate, int chro_same,int gene_parent0[],int gene_parent1[],int gene_baby[]){
    int c0,c1;
    for(int i=0; i<genesize; i++){
        if(i%chro_same == 0){
            c0 = random_int(1);
            c1 = random_int(1);
        }
        gene_baby[i] = mutation_stepwise(gene_parent0[i+c0*genesize], mutation_rate);
        gene_baby[i+genesize] = mutation_stepwise(gene_parent1[i+c1*genesize],mutation_rate);        
    }
}


/* function for mutation (stepwise mutaion) on a locas */
inline int mutation_stepwise(int x,double mutation_rate){
  if(prob(mutation_rate) == 1){
    return x + prob(0.5);
  }
  else return x;
}


/* function for generation of an instance for a new individual */
void birth(int gene_p[2*size_p],int gene_w[2*size_w],int gene_m[2*size_m]){
  int done;
  int count;
  /* if there is dead individual l[i] (from l[0] to l[number_end]), l[i] is recycled for new individual.*/
  done = 0;
  count = 0;
  while((count < number_end+1)&&(done == 0)){
    if(l[count]->alive == 0){   
      l[count]->reset(gene_p, gene_w, gene_m); 
      done = 1;   
      baby_id=count;
    }
    count++;
  } 
  /* if all individual from l[0] to l[number_end] are alive, a new instance l[number_end+k] is
     generated. */
  if(done == 0){
    number_end_buf++;
    l[number_end_buf] = new Ind(gene_p, gene_w, gene_m); 
    baby_id=number_end_buf;
  }
}

//_/_/_/_/_/_/_/_/_/_/_/ functions for setting initial condition  _/_/_/_/_/_/_/_/_/_/_/_/

/* resource distribution for graphics */
void set_R(void){
  double x;
  for(int k=0; k<Z_division; k++){
      x = delta_z*k;
/* combination of three Gausian functions */
      //R[k] = Ra*f_Norm(x, ma, sda)+Rb*f_Norm(x, mb, sdb)+Rc*f_Norm(x, mc, sdc);
      R[k] = R0*f_Norm(x, mR, sdR);
  }  
}

/* initial monomorphic population */
void set_population(int x, int z, int w){
int gene_p0[2*size_p];
int gene_w0[2*size_w];
int gene_m0[2*size_m];

 for(int i=0; i<2*size_p; i++){gene_p0[i] = 0;}
 for(int i=0; i<2*size_w; i++){gene_w0[i] = 0;}
 for(int i=0; i<x; i++){
     gene_p0[i] = 1;
     gene_p0[i+size_p] = 1;
 }
 for(int i=0; i<z; i++){
     gene_w0[i] = 1;
     gene_w0[i+size_w] = 1;
 }

 for(int i=0; i<2*size_m; i++){
     gene_m0[i] = w;   
 }
 
/* Here, setting one ancester individual */
 number_end_buf = number_end;
 birth(gene_p0, gene_w0,gene_m0);
 number_end = number_end_buf;
}




//_/_/_/_/_/_/_/_/ functions for analysis on resulting dynamics _/_/_/_/_/



double species_sx[spe_number];
double species_sy[spe_number];
double species_sx2[spe_number];
double species_sy2[spe_number];
double species_size[spe_number];
int cer[2*size_p+1][amp_size_w*(2*size_w+1)];
int a[4][2] = {{-1,0},{0,1},{1,0},{0,-1}};


/* calculation of phenotype distribution */
void number_count(void){
  number = 0;
  gp=gw=vp=0;
 for(int j=0; j<2*size_p+1; j++){
   for(int k=0; k<amp_size_w*(2*size_w+1); k++){
   phenotype[j][k] = 0;
   }
 }
 maxw=0.0;
  for(int i=0; i<= number_end; i++){
    if(l[i]->alive == 1){
        gp+=l[i]->value_p;
        vp+=l[i]->value_p*l[i]->value_p;
        gw+=l[i]->value_w;
        if(maxw<l[i]->value_w)maxw=l[i]->value_w;
        number++;


      for(int j=0; j<2*size_p+1; j++){
        for (int k=0; k<amp_size_w*(2*size_w+1); k++){
	  if (((j-0.5)*delta_p<= l[i]->value_p)
	      &&(l[i]->value_p < (j+0.5)*delta_p)

              &&((k-0.5)*delta_w <=l[i]->value_w)
	      &&(l[i]->value_w < (k+0.5)*delta_w)){
	    phenotype[j][k]++;
	    //  cerr<<k<<":"<<l[i]->value_w<<" "<<endl;
	  }   
	
	}
      } 
    }
  }
  gp/=double(number);
  gw/=double(number);
  vp/=double(number);
  vp=vp-gp*gp;

}



void cer_wave(int i, int j, int nth,int res){
  int x,y;  
  double gx,gy;
  species_size[nth-1]+=phenotype[i][j];
  species_sx[nth-1]+=(i*delta_p-0.5)*phenotype[i][j];
  species_sy[nth-1]+=j*delta_w*phenotype[i][j];
  species_sx2[nth-1]+=(i*delta_p-0.5)*(i*delta_p-0.5)*phenotype[i][j];
  species_sy2[nth-1]+=j*j*delta_w*delta_w*phenotype[i][j];
  
  if(phenotype[i][j]>wave_edge){
cer[i][j] = nth;
res=resol;
  }
  for(int k=0; k<4; k++){
      x = i+a[k][0];
      y = j+a[k][1];
      if((-1<x)&&( x<2*size_p+1)&&
         (-1<y)&&( y<amp_size_w*(2*size_w+1))&&
         (cer[x][y]==0)){
          if(phenotype[x][y]>wave_edge)cer_wave(x,y,nth,res);
          else{ 
              if(res>0)cer_wave(x,y,nth,res-1);
          }
      }
  }
}


 


void count_species_for_util(void){
  double biomas =0;
  void cer_wave(int i, int j, int nth,int res);
  nth = 1;


for(int i=0; i<spe_number; i++){
species_size[i]=0;
  species_sx[i]=0;
  species_sy[i]=0;  
}
for(int i=0; i<2*size_p+1; i++){
  for(int j=0; j<amp_size_w*(2*size_w+1); j++){
	  cer[i][j]=0;
	}
}
  for(int i=0; i<2*size_p+1; i++){
    for(int j=0; j<amp_size_w*(2*size_w+1); j++){
	  if((phenotype[i][j]>wave_edge)&&(cer[i][j]==0)){
	    cer_wave(i,j,nth,resol);
nth++;	  
}
	  
	  
	}
  }

for(int i=0; i<nth-1; i++){
  species_sx[i] = species_sx[i]/species_size[i];
  species_sy[i] = species_sy[i]/species_size[i];


  species_sx2[i] = species_sx2[i]/species_size[i] - species_sx[i]*species_sx[i];
  species_sy2[i] = species_sy2[i]/species_size[i] - species_sy[i]*species_sy[i];

  biomas +=species_size[i];

}
}
	    	 



void output_util(void){
    cout<<nth-1<<endl;
    for(int i=0; i<nth-1; i++){
        cout<<species_sx[i]<<" "
            << species_sy[i] <<" "
            << species_size[i] <<" "
            <<species_sx2[i]<<" "
            <<species_sy2[i]<<endl;
    }cout<<endl;
}


	  
void output_phenotype_p(void){
    double buf;
   for(int k=1; k<size_p*2+1; k++){
       buf=0;
       for(int j=0; j<amp_size_w*(size_w*2+1); j++)buf+=phenotype[k][j];
       cout<<buf<<" ";   
}
cout<<endl;   

}

void output_phenotype_exist(void){
  int pcount=0;
  double xx0,yy0;
    for(int i=0; i<2*size_p+1; i++){
      for(int j=0; j<amp_size_w*(2*size_w+1); j++){
	if(phenotype[i][j]>0)pcount++;
      }
    }

    cout<<timen<<" "<<pcount<<" "<<2*size_p+1<<" "<<amp_size_w*(2*size_w+1)<<endl;

   for(int i=0; i<2*size_p+1; i++){
      for(int j=0; j<amp_size_w*(2*size_w+1); j++){
	xx0=(i*delta_p-0.5)*phenotype[i][j];
	yy0=j*delta_w*phenotype[i][j];
	if(phenotype[i][j]>0)cout<<xx0<<" "<<yy0<<" "<<phenotype[i][j]<<" "<<i<<" "<<j<<endl;
      }
    }
   cout<<endl;

}


void set_initial_condition(void){
  offset_w=0.0;
    set_R();
    srand(seed);
    number_end = -1;
    /* 10 individuals are introduced for the initial population */
    for(int i=0; i<10; i++)set_population(first_p, first_w,0);
    step = 0;
    number_count();
    for(int i=0;i<number_end+1;i++){
        if(l[i]->alive==1){
            for(int j=0;j<number_end+1;j++){
                if(l[j]->alive==1)l[i]->death+=com(l[i]->value_p,l[j]->value_p)/K(l[i]->value_p);
            }
        }
    }

   for(int k=0; k<size_p*2+1; k++){
     comm[k]=com(delta_p*k,0);
     KK[k]=K(delta_p*k);
   }
}



/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/
/*_/_/_/_/  Functions for graphics and  main loop with X11  _/_/_/_/_/*/
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

#include <X11/Xlib.h>
#include <X11/Xutil.h>



double meshx[g_size_x+1];
double meshy[g_size_y+1];

double deltax=1.0/double(g_size_x);
double deltay=1.0/double(g_size_y);

int g_step=0;
int view_interval = 30;
Display *display;
GC      gc;

XColor iro;
Colormap  cmap;
Window root;
XEvent   event;  

#define MAX_WINDOW 4
Window window[MAX_WINDOW];
Pixmap pix[MAX_WINDOW];
int win_width[MAX_WINDOW],win_height[MAX_WINDOW];
double win_vport[MAX_WINDOW];
bool window_on[MAX_WINDOW]={0};
void (*disp_func[MAX_WINDOW])(void);

int car_id;
bool g_Flag_stop = 0;
bool g_Flag_surface0 = 1;
bool g_Flag_surface1 = 0;
bool g_Flag_water = 0;
bool g_Flag_blend = 0;
bool g_Flag_smooth=1;
bool g_Flag_ortho=0;
bool g_Flag_halt=0;
double g_water_level=1.0;
void g_init(int argc, char *argv[]){
    void set_color(void);
    display = XOpenDisplay(NULL);
    root    = RootWindow(display,0);     
    gc      = XCreateGC(display,root,0,0);
    cmap = DefaultColormap(display, 0);
    set_color();

    for(int i=0; i<g_size_x;i++)meshx[i]=(i/double(g_size_x));
for(int i=0; i<g_size_y;i++)meshy[i]=(i/double(g_size_x));





} 



void g_window2d(int winid, void (*fp)(void),int x, int y, int w, int h){
    win_width[winid]=w;
    win_height[winid]=h;
    window[winid]  = XCreateSimpleWindow(display,root,x,y,w,h,2,0,1);
                                         
    pix[winid]=XCreatePixmap(display,root,w,h,DefaultDepth(display,0));
    XSelectInput(display,window[winid],ButtonPressMask | ExposureMask | KeyPressMask);
    XMapWindow(display,window[winid]);      
    XMoveWindow( display, window[winid], x, y );
    XFlush(display);
    car_id=winid;
    disp_func[winid]=fp;
    window_on[winid]=1;
}


    
unsigned long white,black, green, orange, blue, navy, red, grey;
const int col_grad_number = 30;
unsigned long col_grad_bw[col_grad_number];
unsigned long col_grad_gr[col_grad_number];

unsigned long getcolor(const char* color_name )
     {
       XColor near_color, true_color;
       XAllocNamedColor(display, cmap, color_name, &near_color, &true_color );
       return(near_color.pixel );
     };


unsigned long getrgb(int r, int g, int b)
{
XColor col;    
    col.green = g;
    col.red   = r;
    col.blue  = b;
    XAllocColor( display, cmap, &col );
    return col.pixel;
}


    
void set_color(void){
   white = getcolor("white");
   black = getcolor("black");
   green = getcolor("green");
   orange =  getcolor("orange");
   blue = getcolor("deep sky blue");
   navy = getcolor("navy");
   red = getcolor("red");
   grey = getcolor("dark slate gray");
   for(int i=0; i<col_grad_number; i++){
     col_grad_bw[i] = getrgb(i*65500/col_grad_number,i*65500/col_grad_number,i*65500/col_grad_number);
   }   
   for(int i=0; i<col_grad_number; i++){
     col_grad_gr[i] = getrgb(i*65500/col_grad_number,int(65500*(1-i/double(col_grad_number))),0);
    
   }
   

}

void g_copy(void){
    XCopyArea(display,pix[car_id],window[car_id],
              gc,0,0,win_width[car_id],win_height[car_id],0,0);
    XFlush(display);
};


void g_clear(unsigned long color){
    XSetForeground( display, gc,color);
    XFillRectangle( display, pix[car_id], gc, 0, 0, win_width[car_id], win_height[car_id]);
};

void g_line(double x0,double y0,double x1, double y1){
    
    XDrawLine(display,pix[car_id],gc,int(x0*win_width[car_id]),int(win_height[car_id]-win_height[car_id]*y0),
int(x1*win_width[car_id]),int(win_height[car_id]-win_height[car_id]*y1));
   
}

void g_point(double x0,double y0){
    
  XDrawPoint(display,pix[car_id],gc,int(x0*win_width[car_id]),int(win_height[car_id]-win_height[car_id]*y0));
   
}

void g_line(double x[], double y[],int size,double (*filter)(double y))
{
    
    for (int i = 0; i < size-1; i++) XDrawLine(display,pix[car_id],gc,
int(x[i]*win_width[car_id]),int(win_height[car_id]-win_height[car_id]*filter(y[i])),
int(x[i+1]*win_width[car_id]),int(win_height[car_id]-win_height[car_id]*filter(y[i+1])));
       
}

void g_line_width(int w)
{
    XSetLineAttributes(display,gc,w,LineSolid,CapRound,JoinRound);     
}

void g_color(unsigned long color){
    XSetForeground( display, gc, color);
};


void g_fills(bool (*grad)(int, int)){
    
    for(int i=0; i<g_size_x;i++){
	for(int j=0; j<g_size_y;j++){
if(grad(i,j))XFillRectangle( display, pix[car_id], gc,
int(win_width[car_id]*i/g_size_x), 
win_height[car_id]-int(win_height[car_id]*(j+1)/g_size_y),
int(win_width[car_id]*(i+1)/g_size_x)-int(win_width[car_id]*i/g_size_x),
int(win_height[car_id]*(j+1)/g_size_y)-int(win_height[car_id]*(j)/g_size_y));


	}
    }
}


void g_view_display(void){
    for(int i=0;i<MAX_WINDOW;i++){
        if(window_on[i]){
            car_id=i;
            disp_func[i]();
            
        }
  }
}


    
void g_mouse_2d(int button, int state, int x, int y);


/* Main loop */
void g_main_loop(void){
    
  void keyboard(unsigned char c, int x, int y);
  KeySym key;
  char string[10];
  void g_action(void);
  void g_put_data(void);
  while(!g_Flag_halt){
    
    if(XCheckMaskEvent(display,ButtonPressMask | ExposureMask | KeyPressMask ,&event)){
        switch(event.type)  {
        case Expose:
            
            g_view_display();
            break;
            
        case ButtonPress :        
            g_mouse_2d(event.xbutton.button, event.xbutton.state, event.xbutton.x,event.xbutton.y);
            break;
        case KeyPress:
            XLookupString((XKeyEvent *)&event,string,sizeof(string),&key,NULL);
            keyboard(string[0],0,0);
            break;
        }
        
    }       
    if(!g_Flag_stop){

        g_action();
        g_step++;
        if(g_step%view_interval==0){
            g_put_data();
            g_view_display();
        }
    }
    
    
}
}


/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/
/*_/_/_/   Functions for simulation controll and graphics _/_/_/_/_/_/*/
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

double g_buf[g_size_x];
int g_Flag_3d=1;
int g_Flag_2d=0;

bool mask[g_size_x][g_size_y];
double g_phe[g_size_x][g_size_y];

void g_mouse_2d(int button, int state, int x, int y){
    if (button == 1)  {  
      g_Flag_stop = (g_Flag_stop+1)%2;
      if(g_Flag_stop == 1)cerr<<"Stop!"<<endl;    
    }
    return;
}


void keyboard(unsigned char c, int x, int y) {
}


void g_action(void){
    for(int i=0;i<count_interval;i++)Event();
    number_count();
    
    step++;
    if((Flag_out_species)&&(step%phenotype_interval==0)){
        count_species_for_util();
        cout<<timen<<" "<<nth-1<<endl;
        for(int i=0; i<nth-1;i++){	 
	  cout<<species_sx[i]<<" "<<species_sy[i]+offset_w<<" "<<species_size[i]<<" "<<species_sx2[i]<<" "<<species_sy2[i]<<endl;
        }
	
    }
    
    cerr<<timen<<" "<<nth-1<<" "<<number<<" "<<"male:"<<n_par[0]<<" female:"<<n_par[1]<<endl;    
    
    if((Flag_out_phenotype==1)&&(step%phenotype_interval==0)){
            output_phenotype_exist();
    }
    //        if(timen>10)g_Flag_halt=1;    
}
 
void g_put_data(void){ 
  for(int i=0; i<g_size_x;i++){
      for(int j=0; j<g_size_y;j++){
	if(phenotype[i][j]>0){        
	  mask[i][j]=1;          
	  if(g_Flag_surface1)g_phe[i][j]=phenotype[i][j]*0.02;
	  else g_phe[i][j]=phenotype[i][j]*0.02;
          
	}
	else {
	  mask[i][j]=0;
	  g_phe[i][j]=0.0;
	}
      }
  }
  
}


inline double filter2(double y){
    return y*1.0e-4;
}

inline bool grad1(int i, int j){
    double col;
    int pi;

    if(phenotype[i][j]>0){
        col=phenotype[i][j];       
        g_color(col_grad_gr[min(int(col),col_grad_number-1)]);
        return 1;
    }
    else{
       return 0;
    }
        
}
    

void display1(void){
    g_clear(navy);
    g_fills(grad1);
    g_color(green);
    g_line(meshx,R,g_size_x,filter2);
    g_copy();

}

void display_generation(int argc, char **argv){
    view_interval=1;
    set_initial_condition();
    g_water_level = 0.2;
    g_Flag_blend = 1; 
    g_Flag_water = 0;
    g_Flag_surface0 = 1;    
    g_init(argc,argv);
    g_window2d(window[0],display1,100,0,400,amp_size_w*200);
    g_main_loop();
      
}

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/  Main  _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

main(int argc, char *argv[]){
  int opt;
    while((opt=getopt(argc, argv,
                    "pu")) != EOF){
    switch(opt){
    case 'u':
      Flag_out_species=1;
      break;
    case 'p':
    Flag_out_phenotype=1;
      break;
      
       
    }
    }
  
    display_generation(argc,argv);
}
