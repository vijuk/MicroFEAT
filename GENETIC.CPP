
/*   VARIABLE  DECLARATIONS  */
  
# define NEG  0
# define ZERO 1
# define POS  2

int  initialize();
void generation();
void statistics();
void report();
void initpop();
int  select();
int  mutation();
void crossover();
void part_crossover();
void initial_fitness();
							   10101010100101
typedef struct				   
{
  unsigned int gene : 1;
} chromosome;

typedef struct
{
   chromosome chrome[120];
   int cntrl_set[5],T_op;
   double control_peak[5],control_width[5];
   double error_peak[5],error_width[5];
   double ce_peak[5],ce_width[5];
   double fitness;
} individual;

individual *oldpop,*newpop;    
int popsize,lchrome,maxgen;
int c,mate1,mate2,no_of_rules = 4;
double pcross,pmutation,sumfitness;
double oldavg,newavg,newmax,oldmax,min;

void decode();
double objfunct();

/*  RANDOM   */

void randomize();
double random();
int flip();
int rnd();
void advance_random();
void warmup_random();
int trunc();

double oldrandom[55];
int jrand;

/* GENETIC ALGORITHM FOR STRUCTURAL OPTIMIZATION  */

# include <stdio.h>
# include <math.h>

FILE *f;
void main()
{
   int gen,fit,i,j,repeat;
   double yy;

   f = fopen("data.out","w");
   
   printf("\n\n\n      G A  FOR FUZZY CONTROLLER DESIGN \n\n");
   printf("\n1 ->  first order\n");
   printf("2 ->  second order\n");
   printf("3 ->  first order with delay\n");
   printf("\nSelect the system   --> ");
   
   scanf("%d",&c);
   if(c == 1 || c == 2 || c == 3)
   {
      fit = 0; repeat = 0;
      if(initialize() == 0) 
      {
	  outtextxy(30,240,"\n\n Memory allocation error!!!");
	  getch();
	  return;
      }

      printf("\n\n process going on .......");
      for(gen=1;gen <= maxgen;gen++)
      {
	   for(i=0;i<popsize;i++)
              oldpop[i] = newpop[i];  
           oldavg = newavg;
           oldmax = newmax;
           generation();
	   statistics();
	   if(newavg/oldavg <= 1.0 && newmax <= oldmax) 
           {
	      repeat ++;
	      if(repeat == 2) {  gen ++;  break; }
           }
           else repeat = 0;
      }
      for(i=0;i<popsize;i++)
        if(newpop[i].fitness > newpop[fit].fitness) fit = i;
      fprintf(f,"  Generation is  %d\n",gen);
      fprintf(f,"  Maximum  fit string is\n");
      for(i=0;i<lchrome;i++)
	fprintf(f,"%d",newpop[fit].chrome[i].gene);
      fprintf(f,"\n  Fitness   =  %lf",newpop[fit].fitness);
      fprintf(f,"\n\n");
      printf("\npress a key   ");
      getch();
	 outtextxy(50+50*i,402,"*");
   }
   fclose(f);
   getch();
   closegraph();
}   /*  end of main  */

int initialize()  /* read  parameters */
{
   int i;
   cleardevice();
   gotoxy(2,4);
   printf("\n");
   do
   {
      printf("\n Enter population size   -->  ");
      scanf("%d",&popsize);
   }
   while(popsize <= 0 || popsize/2*2 != popsize);
   do
   {
      printf("\n Enter max. generations   -->  ");
      scanf("%d",&maxgen);
   }
   while(maxgen <= 0);
   lchrome = 106;
   do
   {
      printf("\n Enter crossover probability   -->  ");
      scanf("%lf",&pcross);
   }
   while(pcross>1 || pcross<0); 
   do
   {
      printf("\n Enter mutation  probability   -->  ");
      scanf("%lf",&pmutation);
   }
   while(pmutation>1 || pmutation<0); 
   if((oldpop = (individual *)malloc((popsize+1)*sizeof(individual))) == 0)
       return 0;
   if((newpop = (individual *)malloc((popsize+1)*sizeof(individual))) == 0)
       return 0;
   randomize();
   cleardevice();
   initpop();
   cleardevice(); 
   gotoxy(4,4); 
   return 1;
}

void initpop() /* initialises the population */
{
   int j,k;
   char  n;
   FILE  *fp;
   fp = fopen("fuzzy.dat","r");
   for(k=0;k<lchrome;k++)  
   {
        n = fgetc(fp); 
        if(n == '1') newpop[0].chrome[k].gene = 1;
	else  newpop[0].chrome[k].gene = 0;
   } 
   decode(0,newpop[0].chrome); 
   fclose(fp);
   do
   {
      fflush(stdin);
      printf("\n Whether to consider the initial rule (y/n)  -> ");
      scanf("%c",&n);
      if((n = toupper(n)) == 'Y') j = 1;
      else if(n == 'N') j = 0;
   } while(!(n == 'Y' || n == 'N'));
   for(;j<popsize;j++)  
   {
         for(k=0;k<lchrome;k++)  
         {
	    if(flip(0.5)) newpop[j].chrome[k].gene = 1;
	    else  newpop[j].chrome[k].gene = 0;
         } 
         decode(j,newpop[j].chrome); 
   } 
   cleardevice();
   gotoxy(4,4);
   printf("\n process going on ....");
   for(j=0;j<popsize;j++)  
      newpop[j].fitness = objfunct(j);
   cleardevice();
   gotoxy(4,4);
   statistics();
   initial_fitness();
   statistics();
   printf("\n\n\n\n");
   printf("\n     Initial max. fitness = %lf",newmax);
   printf("\n     Initial min. fitness = %lf",min);
   printf("\n     Initial avg. fitness = %lf",newavg);
   printf("\n     Initial sum. fitness = %lf",sumfitness);
   printf("\n\n   Press <return> to continue ");
   getch();
   return;
}
 
void generation()  /* find the next generation */
{
   int j;								
   for(j=0;j<popsize;j+=2)				
   {
      mate1 = select();
      mate2 = select();
      crossover(j);
      decode(j,newpop[j].chrome);
      decode(j+1,newpop[j+1].chrome);
   }
   for(j=0;j<popsize;j++)
      newpop[j].fitness =  objfunct(j);
   return;
}

int select()   /* roulette wheel selection of chromosomes */
{
   double ran,partsum=0;
   int j=0;
   ran = random()*sumfitness;
   do 
   {
      partsum += oldpop[j].fitness;
      j ++;
   } while(partsum<ran && j<popsize);
   return (j-1);
}

int  mutation(int i,int j)
{
   chromosome c;
   if(i == 1) i = mate1;
   else i = mate2;
   if(flip(pmutation))
      c.gene = ~oldpop[i].chrome[j].gene;
   else  c.gene = oldpop[i].chrome[j].gene;
   return  c.gene;
}

void crossover(int ipop)
{
   int j,start,end;
   for(j=0;j<no_of_rules-1;j++)
   {
      start = j*26+0;  end = j*26+25;
      part_crossover(start,end,ipop);
   }
   start = j*26+0;  end = j*26+27;
   part_crossover(start,end,ipop);
   return;
}

void part_crossover(int start,int end,int ipop)
{
   int j,jcross,k;
   if(flip(pcross)) 
      jcross = rnd(start,end);
   else jcross = end;
   for(j=0;j<=jcross;j++)
   {
      newpop[ipop].chrome[j].gene = mutation(1,j);
      newpop[ipop+1].chrome[j].gene = mutation(2,j);  
   }  
   for(j=jcross+1;j<=end;j++)
   {
      newpop[ipop].chrome[j].gene = mutation(2,j);
      newpop[ipop+1].chrome[j].gene = mutation(1,j);
   }
   return;
}

void statistics()   /*  find average value  */
{
   int j=0;
   sumfitness = newpop[0].fitness;
   newmax = newpop[0].fitness;
   min = newpop[0].fitness;
   for(j=1;j<popsize;j++)
   {
      sumfitness += newpop[j].fitness;
      if(newpop[j].fitness > newmax) newmax = newpop[j].fitness;
      if(newpop[j].fitness < min) min = newpop[j].fitness;
   }
   newavg = sumfitness/(float)popsize;
   return;
}

void report(int gen)    /*  print the chromosome details  */
{
   int i,k=0,j,set;
   cleardevice();
   printf("\n  GENERATION  -->  %d\n\n\n\n",gen);
   for(set=0;set<no_of_rules;set++)
   {
     printf("   RULE SET   %d\n",set+1);
     for(i=0;i<popsize;i++)
     {
        printf("chromosome[%d] = ",i+1);
        for(j=0;j<26;j++)
           printf("%d",newpop[i].chrome[j+set*26].gene);
        printf("   Fitness  = %lf\n\n",newpop[i].fitness);
        if(++k == 7)
        {
          printf("\nPress <return> to continue  ");
	  getch();
	  k = 0;  
	  cleardevice();
	  printf("\n  GENERATION  -->  %d\n\n\n\n",gen);
	}
     }
   }
   printf("\n\nPress <return> to continue  ");
   fflush(stdin);
   getch();
   return;
}

void initial_fitness()
{
  int i;
  double range = newmax-min;
  for(i=0;i<popsize;i++)
    newpop[i].fitness = min + 5.0*newpop[i].fitness/range;
  return;
}
      /*    random.funct  */

  /*   RANDOM GENERATION  */  

void randomize()
{
   double randomseed;
   do
   {
     printf("\n Enter seed random no.(0 .. 1)  -->  ");
     scanf("%lf",&randomseed);
   } while(randomseed<=0 || randomseed>=1);
   warmup_random(randomseed);
   return;
}

double random()  /* returns a value between 0 and 1  */
{
   jrand ++;
   if(jrand > 54)
   {
      jrand = 0; advance_random();    
   }
   return oldrandom[jrand];
}

int flip(double prob)  /* returns 0 or 1 with probability = prob */
{
  if(prob == 1.0) return 1;
  else
  if(random()<=prob)  return 1;
  else return 0;
}

int rnd(int low,int high) /* returns a number between low and high  */
{
   int i;
   if(low >= high) i = low;
   else
   {
      i = trunc(random()*(high-low+1)+low);
      if(i > high) i = high;
   }
   return i;
}

void advance_random()
{
   int j1;
   double newrandom;
   for(j1=0;j1<=23;j1++)
   {
      newrandom = oldrandom[j1]-oldrandom[j1+31];
      if(newrandom < 0.0) newrandom ++;
      oldrandom[j1] = newrandom;
   }
   for(j1=24;j1<=54;j1++)
   {
      newrandom = oldrandom[j1]-oldrandom[j1-24];
      if(newrandom < 0.0) newrandom ++;
      oldrandom[j1] = newrandom;
   }
   return;
}

void warmup_random(double random_seed)
{
   int i,j;
   double  newrandom,prevrandom;
   oldrandom[54] = random_seed;
   newrandom = 1.0e-9;
   prevrandom = random_seed;
   for(j=1;j<=54;j++)
   {
       i = (21*j)%55;
       oldrandom[--i] = newrandom;
       newrandom = prevrandom-newrandom;
       if(newrandom < 0.0) newrandom ++;
       prevrandom = oldrandom[i];
   }
   advance_random();
   advance_random();
   advance_random();
   jrand = 0;
   return;
}

int trunc(double y)
{
   return (int)(y+0.5);
}
