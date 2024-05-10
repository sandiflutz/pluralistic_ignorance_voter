/***************************************************************************************************************************************************** 
 *           -This program was to explore pluralistic ignorance.                                                                                     *
 *           -This program simulates a voter model with concealed opinions (public and private opinions=external and internal opinions).             *  
 *           -The networks used here are: well-mixed with groups and square lattice.                                                                 *
 *           -There are 2 possible opinions: 0 (lack of action, for external states, or lack of intention to act, for internal states) or            *
 *           1 (action, for external states, or intention, for internal states).                                                                     *
 *           -The dynamics consists of imitation, externalization of internal opinions and internalization of external opinions events               *
 *           -Agents have a memory of length T_mem. In the internalization events, they internalize the most frequent external state they chose      *
 *           in their last T_mem oportunities                                                                                                        *    
 *****************************************************************************************************************************************************/


#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mc.h>
#include<assert.h>

/**********defining global constants and macros********************/

#define L                  22  /*square lattice linear size*/
#define N                  (L*L)/*square lattice size*/
#define VIZ                20   /*number of nearest neighbors each agent has (lattice version)*/
#define Gsize              20   /*number of agents in each group (well mixed version)*/
#define Pint               0.01 /*probability of internalizing an external opinion*/
#define Pcp                0.5/*probability of copying the majorities strategy*/
#define GAMMA              0.7 /*exponent for the probability of copying the external opinion 1*/

#define Di1                0.5 /*fraction of internal opinioms that are 1*/

#define REDE               0 /*type of network. 0: well mixed with groups; 1: square lattice with Viz nearest neighbors*/
#define INIT               1/*0: random (uniform); 
			     *1: external ops. start as 0, while internal ops. are randomly chosen*/
#define T_mem              1/*memory size in time steps*/
#define Teq                10000 /*time for the system to reach equilibrium*/
#define Tf                 (Teq+100)
#define SAVE_CONFIG_ID     0  /*0:png
			       *1:eps
			       */
#define SIZE_CORR         
#define SAMPLE             100
/*****************declaring routines order**************************/
/*----initialization routines-----*/
void initialization(void);
void initialState(void);
void setNeighbors(void);
/*---dynamics-----*/
void evolve(int tf);
void dynamics(void);
int randomNeighbor(int x);
/*--measures----*/
void calcDens(void);
void save_densXt(int tf);
void densXpint(int tmax);
void densXpcp(int tmax);
void corrXpcp(int tmax);
void densXgama(int tmax);
void diagram_pcpXpintXd(int tmax);
void diagram_pcpXgamaXd(int tmax);
void save_config(int t);
/*----finalization-----*/
void freeMemory(void);
/*****************Global Variables******************************/
int numsteps, tmem, intern,gsize,*s[2],*s0[2],*list,*mem, *mem0,*group, **neighbor,*count_mem,*s_old[2];
double dens[2],dens2[4],pint,pcp,gama,dens_ext1,dens_ext0,dens_int0,dens_int1,dens_cp1,dens_cp0,dens_mem0,dens_mem1,*corr_e,*corr_i,*corr_ie[4];
FILE *fdensXt,*fdensXpint,*fdensXpcp,*fdensXgama,*fdiagpcpi,*fdiagpgd,*fdiagppd,*fcorrXpcp;
/****************Program's Routines*****************************/

int main(void){

	initialization();

//--------------------calling measure routines--------------------//

#ifdef DENSxT
	int i;
	for(i=0; i<1; ++i){
		save_densXt(Tf);
	}
#endif
#ifdef DENSxPint
	int i;
	for(i=0; i<SAMPLE; ++i){
		densXpint(Tf);
	}
#endif
#ifdef DENSxPcp
	int i;

	for(i=0; i<SAMPLE; ++i){
		densXpcp(Tf);
	}
#endif
#ifdef CORRxPcp
	int i;

	for(i=0; i<SAMPLE; ++i){
		corrXpcp(Tf);
	}
#endif
#ifdef DENSxGAMA
	int i;

	for(i=0; i<SAMPLE; ++i){
		densXgama(Tf);
	}
#endif
#ifdef DIAG_PcpXPintXD
	int i;
	for(i=0; i<SAMPLE; ++i){
		diagram_pcpXpintXd(Tf);
	}
#endif
#ifdef DIAG_PcpXGAMAXD
	int i;
	for(i=0; i<SAMPLE; ++i){
		diagram_pcpXgamaXd(Tf);
	}
#endif
#ifdef SAVE_CONFIG
	evolve(2000);
#endif
//--------free allocated memory and closes files that are still open--------/

	freeMemory();
	return 0;
}
/*******************************************************
 *   initializes most important features of the system *
 *   (lattice and neighboring structure,               *
 *   global variables, etc.)                           *
 ******************************************************/
void initialization(void){
	int i;
	unsigned int seed;
	seed=start_randomic(0);


	pint=Pint;
	pcp=Pcp;
	gama=GAMMA;
	tmem=T_mem;
	gsize=Gsize;

	for(i=0; i<2; ++i){//s[0][j]: internal op. of j; s[1][j]: external op. of j;
		s[i]=(int *)malloc(sizeof(int)*N);
	}

	count_mem=(int *)malloc(sizeof(int)*N);

	mem=(int *)malloc(sizeof(int)*N);
	mem0=(int *)malloc(sizeof(int)*N);

	neighbor=(int **)malloc(sizeof(int *)*VIZ);
	for(i=0; i<VIZ; ++i){
		neighbor[i]=(int *)malloc(sizeof(int)*N);
	}
#if  (REDE==0)
	list=(int *)malloc(sizeof(int)*N);//list of positions of the available players to choose to be part of a group
	group=(int *)malloc(sizeof(int)*gsize);
#endif
#if (REDE==1)
	setNeighbors();
#endif
	for(i=0; i<N; ++i){
		mem0[i]=(int)(FRANDOM*2.);
	}
	initialState();

	return;
}
/*****************************************************
 *   set initial conditions                          *
 *****************************************************/
void initialState(void){
	int i,j;
	double center,dist,r;

#if  (INIT==0)//random
	for(i=0; i<N; ++i){
		for(j=0; j<2; ++j){//j=0 is the internal op., j=1 is external op.
			s[j][i]=(int)(FRANDOM*2.);
		}
	}
#endif
#if  (INIT==1)//
	for(i=0; i<N; ++i){
		s[1][i]=0;//external opinions start as 0
		if(FRANDOM<Di1){
			s[0][i]=1;
		}else{
			s[0][i]=0;
		}
	}
#endif
	return;
}
/***********************************************************************
 *    set matrix of neighbors for the square lattice version           *
 ***********************************************************************/
void setNeighbors(void){
	int i,k,base;

	k=0;
	while(k<VIZ){
		for(i=0; i<N; ++i){
			base=(i/L)*L;
			switch(k){
				//4 neighbors
				case 0: neighbor[k][i]= base + (i-1+L)%L;//left
					break;
				case 1: neighbor[k][i]= (i-L+N)%N;//up
					break;
				case 2: neighbor[k][i]= base + (i+1)%L;//right
					break;
				case 3: neighbor[k][i]= (i+L)%N;//down
					break;
 			#if (VIZ>=8)
				case 4: //right-up
					neighbor[k][i]= base + (i+1)%L;//(right)
					neighbor[k][i]=(neighbor[k][i]-L+N)%N;//(up)
					break;
				case 5: //right-down
					neighbor[k][i]= base + (i+1)%L;//(right)
					neighbor[k][i]= (neighbor[k][i]+L)%N;//(down)
					break;
				case 6: //left-up
					neighbor[k][i]= base + (i-1+L)%L;//(left)
					neighbor[k][i]= (neighbor[k][i]-L+N)%N;//(up)
					break;
				case 7: //left-down
					neighbor[k][i]= base + (i-1+L)%L;//(left) 
					neighbor[k][i]= (neighbor[k][i]+L)%N;//(down)
					break;
 			#endif
			#if (VIZ>=12)
				case 8: //left-left
				       	neighbor[k][i]= base + (i-1+L)%L;//(left)
					neighbor[k][i]= (neighbor[k][i]-1+L)%L;//(left)
					break;
				case 9: //up-up
					neighbor[k][i]= base + (i-L+N)%N;//(up)
					neighbor[k][i]= (neighbor[k][i]-L+N)%N;//(up)
					break;
				case 10://right-right 
					neighbor[k][i]=base + (i+1)%L;//(right)
					neighbor[k][i]=(neighbor[k][i]+1)%L;//(right)
					break;
				case 11://down-down 
					neighbor[k][i]= base +(i+L)%N;//(down)
					neighbor[k][i]= (neighbor[k][i]+L)%N;//(down)
					break;
			#endif
			#if (VIZ>=20)
				case 12: //left-left-up 
					neighbor[k][i]= base + (i-1+L)%L;//(left)
					neighbor[k][i]= (neighbor[k][i]-1+L)%L;//(left)
					neighbor[k][i]=(neighbor[k][i]-L+N)%N;//(up)
					break;
				case 13://left-left-down 
					neighbor[k][i]= base + (i-1+L)%L;//(left) 
					neighbor[k][i]= (neighbor[k][i]-1+L)%L;//(left)
					neighbor[k][i]=(neighbor[k][i]+L)%N;//(down)
					break;
				case 14://right-right-up 
					neighbor[k][i]= base + (i+1)%L;//(right) 
					neighbor[k][i]=(neighbor[k][i]+1)%L;//(right)
					neighbor[k][i]=(neighbor[k][i]-L+N)%N;//(up) 
					break;
				case 15://right-right-down 
					neighbor[k][i]= base + (i+1)%L;//(right) 
					neighbor[k][i]=(neighbor[k][i]+1)%L;//(right)
					neighbor[k][i]=(neighbor[k][i]+L)%N;//(down)
					break;
				case 16://up-up-left 
					neighbor[k][i]= base + (i-L+N)%N;//(up) 
					neighbor[k][i]= (neighbor[k][i]-L+N)%N;//(up)
					neighbor[k][i]= (neighbor[k][i]-1+L)%L;//(left)
					break;
				case 17://up-up-right 
					neighbor[k][i]= base + (i-L+N)%N;//(up) 
					neighbor[k][i]= (neighbor[k][i]-L+N)%N;//(up)
					neighbor[k][i]=(neighbor[k][i]+1)%L;//(right)
					break;
				case 18://down-down-left 
					neighbor[k][i]= base + (i+L)%N;//(down) 
					neighbor[k][i]= (neighbor[k][i]+L)%N;//(down)
					neighbor[k][i]= (neighbor[k][i]-1+L)%L;//(left)
					break;
				case 19://down-down-right 
					neighbor[k][i]= base + (i+L)%N;//(down) 
					neighbor[k][i]= (neighbor[k][i]+L)%N;//(down)
					neighbor[k][i]=(neighbor[k][i]+1)%L;//(right)
					break;
			#endif
			}
		}
		++k;
	}
	return;
}
/*******************************************
*   evolves the system in time             *
*******************************************/
void evolve(int tf){
	int i,j,id,ne,ni,soma;
	int save_this,nf,num_files,interv_grav;//necessary for save_config routine


#ifdef SAVE_CONFIG//save spatial configuration (for lattice version)
	nf=0;
	#if (SAVE_CONFIG_ID==0)//multiple files (png) for creating an animation
	interv_grav=1;
	save_this=0;
	num_files=(int)tf/interv_grav;
	#else //files for creating eps snapshots
	interv_grav=1;
	save_this=tf;
	num_files=1;
	#endif
#endif


	/*set initial memory*/
	for(i=0; i<N; ++i){
		mem[i]=mem0[i];//random memory
		count_mem[i]=0;
	}
	/*-------------------*/

/*-----------time loop------------------*/
	numsteps=0;
	while(numsteps<=tf){
		dynamics();
		calcDens();
#ifdef DENSxT//store densities as functions of time
		fprintf(fdensXt,"%d %f %f %f %f\n",numsteps,dens2[0],dens2[1],dens2[2],dens2[3]);
		fflush(fdensXt);
		printf("%d %f %f %f %f\n",numsteps,dens2[0],dens2[1],dens2[2],dens2[3]);
#endif
#ifdef SAVE_CONFIG//store spatial configurations
		if((numsteps==save_this)&&(nf<=num_files)){
			save_config(save_this);
			save_this += interv_grav;
			++nf;
			printf("%d\n",numsteps);
		}
#endif
		++numsteps;
	}

	return;
}
/*******************************************
 *   imitation, externalization and        *
 *   internalization dynamics              *
 ******************************************/
void dynamics(void){
	int i,j,k,id,x,viz,num_gr,sum_maioria,n_list,temp;
	double prob_cp1,tmeio;



#if  (REDE==0)//well mixed with groups 
		
	for(j=0; j<N; ++j){
		for(i=0; i<N; ++i){
			list[i]=i;//initialization of the position list of available agents in the system (which will be randomly choose to be part of a group)
		}
		n_list=N;
		for(i=0; i<gsize; ++i){//randomly chooses available agents (=that are not part of a group) in the system to be part of a group
			x=(int)((double)FRANDOM*n_list);//lembrar de testar se x já foi escolhido
			group[i]=list[x];
			list[x]=list[n_list-1];
			--n_list;	
		}
//----------GROUP-DYNAMICS------------------------------------------------------------------------------------------//
		sum_maioria=0;
		for(i=0; i<gsize; ++i){
//-----------imitation and externalization events for the random groups simulation--------------//
			prob_cp1=(double)sum_maioria/gsize;
			prob_cp1=pow(prob_cp1,gama);
			if(FRANDOM<pcp){//imitates
				if(FRANDOM<prob_cp1){//imitates external state 1
					s[0][group[i]]=1;
					++sum_maioria;
				}else{//imitates external state 0
					s[0][group[i]]=0;	
				}
			}else{//externalizes internal opinion
				s[0][group[i]]=s[1][group[i]];
				sum_maioria+=s[0][group[i]];
			}


//-----------memory update for the random groups simulation------------------------//
			if(s[0][group[i]]==0){
				--count_mem[group[i]];
				if(count_mem[group[i]]<0)count_mem[group[i]]=0;//count_mem cannot be negative
			}else{
				++count_mem[group[i]];
				if(count_mem[group[i]]>T_mem)count_mem[group[i]]=T_mem;//count_mem cannot be greater than T_mem
			}

			tmeio=T_mem/2.;
			if(count_mem[group[i]]>tmeio){
				mem[group[i]]=1;
			}else if(count_mem[group[i]]<tmeio){
				mem[group[i]]=0;
			}else{//this happens only when T_mem is even, otherwhise tmeio won't be an integer
				mem[group[i]]=(int)(FRANDOM*2.);
			}
							    
//-----------------internalization for the random groups simulation------------------//
			if(FRANDOM<pint){//memory is internalized with prob pint
				s[1][group[i]]=mem[group[i]];
		
			}
		}

	}

#endif

//----------SQUARE-LATTICE-DYNAMICS------------------------------------------------------------------------------------------//
#if  (REDE==1)//square lattice with 4 nearest neighbors
	for(i=0; i<N; ++i){
		x=(int)(FRANDOM*N);//x is the position of agent in the center of the group
		sum_maioria=0;
//-----------imitation and externalization events for the square lattice simulation--------------//
		if(FRANDOM<pcp){//agent x imitates
			prob_cp1=(double)sum_maioria/(VIZ+1);//prob_cp1 is the probability of imitating state 1
			prob_cp1=pow(prob_cp1,gama);
			if(FRANDOM<prob_cp1){//agent x imitates state 1
				s[1][x]=1;
				++sum_maioria;
			}else{//agent x imitates state 0
			        s[1][x]=0;
			}
		}else{//agent x externalizes internal opinion
			s[1][x]=s[0][x];
			sum_maioria+=s[1][x];
		}


//-----------memory update for the square lattice simulation------------------------//
		if(s[1][x]==0){
			--count_mem[x];
			if(count_mem[x]<0)count_mem[x]=0;//count_mem cannot be negative
		}else{
			++count_mem[x];
			if(count_mem[x]>T_mem)count_mem[x]=T_mem;//count_mem cannot be grater than T_mem
		}
		tmeio=T_mem/2.;
		if(count_mem[x]>tmeio){
			mem[x]=1;
		}else if(count_mem[x]<tmeio){
			mem[x]=0;
		}else{//this happens only when T_mem is even, otherwhise tmeio won't be an integer
			mem[x]=(int)(FRANDOM*2.);
		}
//-----------------internalization for the random groups simulation------------------//
		if(FRANDOM<pint){
			s[0][x]=mem[x];//internalizes memory
		}
//-----------imitation, externalization and internalization events for the neighbors of x in the square lattice simulation--------------//
		for(j=0; j<VIZ; ++j){//neighbors of agent in x
			viz=neighbor[j][x];
//-----------imitation and externalization-----------------------------/
			if(FRANDOM<pcp){//imitation
				prob_cp1=(double)sum_maioria/(VIZ+1);
				prob_cp1=pow(prob_cp1,gama);
				if(FRANDOM<prob_cp1){//imitates state 1
					s[1][viz]=1;
					++sum_maioria;
				}else{//imitates state 0
					s[1][viz]=0;
				}
			}else{//externalization
				s[1][viz]=s[0][viz];
				sum_maioria+=s[1][viz];
			}

//-----------memory update for the neighbors of x in the square lattice simulation------------------------//
			if(s[1][viz]==0){
				--count_mem[viz];
				if(count_mem[viz]<0)count_mem[viz]=0;//count_mem nao pode ser negativo
			}else{
				++count_mem[viz];
				if(count_mem[viz]>T_mem)count_mem[viz]=T_mem;//count_mem nao pode ser maior que T_mem
			}
			tmeio=T_mem/2.;
			if(count_mem[viz]>tmeio){
				mem[viz]=1;
			}else if(count_mem[viz]<tmeio){
				mem[viz]=0;
			}else{//happens only for an even T_mem, otherwise tmeio wont be and interger 
				mem[viz]=(int)(FRANDOM*2.);
			}
//-----------internalization for the neighbors of x in the square lattice simulation------------------------//
			if(FRANDOM<pint){
				s[0][viz]=mem[viz];
			}
		}

	}
#endif
		return;
}
/******************************************
 * selects a random neighbor of agent in  *
 * position x                             *
 *****************************************/
int randomNeighbor(int x){
	int y,id;

	id=(int)(FRANDOM*(double)VIZ);
	y=neighbor[id][x];

	return y;
}
/*********************************************
 *  calculates density of opnion 1           *
 ********************************************/
void calcDens(void){
	int i,j,id;

	dens[0]=0.; 
	dens[1]=0.;
	memset(dens2,0.,sizeof(float)*4);
	for(i=0; i < N; ++i){
		id=2*s[0][i]+s[1][i];//id: 00->0, 01->1, 10->2, 11->3
		dens2[id]+=1.;
		dens[0]+=(double)s[0][i];
		dens[1]+=(double)s[1][i];
	}
	dens[0]/=(double)N;
	dens[1]/=(double)N;
	for(i=0; i<4; ++i){
		dens2[i]/=(double)N;
	}
	return;
}
/*********************************************
 *    Save the density of both internal and  *
 *    external opinions as functions of time *
 *********************************************/
void save_densXt(int tmax){
	int i,ok;
	char name[100];
	unsigned long id;

	//creating multiple files for calculating averages outside of the program 
	id = (unsigned long)time(NULL);

	ok=0;
	while(ok==0){
	#if (REDE==0)
		sprintf(name,"densXt_L%d_Gsize%d_G%0.2f_Pcp%0.2f_Pext%0.2f_Pint%0.2f_Di1%0.2f_Tmem%d_%ld.dat",L,gsize,gama,pcp,pext,pint,Di1,T_mem,id);
	#else
		sprintf(name,"densXt_L%d_Viz%d_G%0.2f_Pcp%0.2f_Pext%0.2f_Pint%0.2f_Di1%0.2f_Tmem%d_%ld.dat",L,VIZ,gama,pcp,pext,pint,Di1,T_mem,id);
	#endif
		fdensXt = fopen(name,"r");
		if(fdensXt!=NULL){//testing if a file with the above name already exists 
			++id; 
			fclose(fdensXt);
		}else{
			ok=1;
		}
	}

	#if  (REDE==0)
	sprintf(name,"densXt_L%d_Gsize%d_G%0.2f_Pcp%0.2f_Pext%0.2f_Pint%0.2f_CI%d_Di1%0.2f_Tmem%d_%ld.dat",L,gsize,gama,pcp,pext,pint,INIT,Di1,T_mem,id);
	#else
	sprintf(name,"densXt_L%d_Viz%d_G%0.2f_Pcp%0.2f_Pext%0.2f_Pint%0.2f_CI%d_Di1%0.2f_Tmem%d_%ld.dat",L,VIZ,gama,pcp,pext,pint,INIT,Di1,T_mem,id);
	#endif
	fdensXt = fopen(name,"w");


	fprintf(fdensXt,"#1:time 2:dens[0] 3:dens[1]\n");
	///////////////////////////

	initialState();

	evolve(tmax);

	fclose(fdensXt);
	return;
}
/***********************************************
 *      densities as functions of the          *
 *      internalization probability            *
 **********************************************/
void densXpint(int tmax){
	int i,ok;
	double dp,pmax;
	unsigned long id;
	char name[100];

	id = (unsigned long)time(NULL);

	ok=0;
	while(ok==0){
	#if  (REDE==0)
		sprintf(name,"densXpint_L%d_Gsize%d_G%0.2f_Pcp%0.2f_Pext%0.2f_%ld.dat",L,gsize,gama,pcp,pext,id);
	#else
		sprintf(name,"densXpint_L%d_VIZ%d_G%0.2f_Pcp%0.2f_Pext%0.2f_%ld.dat",L,VIZ,gama,pcp,pext,id);
	#endif
		fdensXpint = fopen(name,"r");
		if(fdensXpint!=NULL){
			++id;
			fclose(fdensXpint);
		}else{
			ok=1;
		}
	}

	#if  (REDE==0)
	sprintf(name,"densXpint_L%d_Gsize%d_G%0.2f_Pcp%0.2f_Pext%0.2f_%ld.dat",L,gsize,gama,pcp,pext,id);
	#else
	sprintf(name,"densXpint_L%d_VIZ%d_G%0.2f_Pcp%0.2f_Pext%0.2f_%ld.dat",L,VIZ,gama,pcp,pext,id);
	#endif
	fdensXpint = fopen(name,"w");

	fprintf(fdensXpint,"#1:pint 2:dens00 3:dens01 4:dens10 5:dens11\n");
	dp=0.01;
	pint=dp;
	pmax=0.9;
	while(pint<=pmax){
		initialState();
		evolve(tmax);

		fprintf(fdensXpint,"%f ",pint);
		fflush(fdensXpint);
		 printf("%f ",pint);
		for(i=0; i<4; ++i){
			fprintf(fdensXpint,"%f ",dens2[i]);
			fflush(fdensXpint);
			printf("%f ",dens2[i]);
		}
		fprintf(fdensXpint,"\n");
		fflush(fdensXpint);
		printf("\n");

		pint+=dp;
	}


	fclose(fdensXpint);
	return;
}
/***********************************************
 *      densities as functions of the          *
 *      imitation probability                  *
 **********************************************/
void densXpcp(int tmax){
	int i,ok,num[4];
	double dp,pmax;
	unsigned long id;
	char name[100];

	id = (unsigned long)time(NULL);

	ok=0;
	while(ok==0){
	#if  (REDE==0)
		sprintf(name,"densXpcp_L%d_Gsize%d_G%0.2f_Pext%0.2f_Pint%0.2f_CI%d_Tmem%d_%ld.dat",L,gsize,gama,pext,pint,INIT,T_mem,id);
	#else
		sprintf(name,"densXpcp_L%d_VIZ%d_G%0.2f_Pext%0.2f_Pint%0.2f_CI%d_Tmem%d_%ld.dat",L,VIZ,gama,pext,pint,INIT,T_mem,id);
	#endif
		fdensXpcp = fopen(name,"r");
		if(fdensXpcp!=NULL){
			++id;
			fclose(fdensXpcp);
		}else{
			ok=1;
		}
	}

	#if  (REDE==0)
	sprintf(name,"densXpcp_L%d_Gsize%d_G%0.2f_Pext%0.2f_Pint%0.2f_CI%d_Tmem%d_%ld.dat",L,gsize,gama,pext,pint,INIT,T_mem,id);
	#else
	sprintf(name,"densXpcp_L%d_VIZ%d_G%0.2f_Pext%0.2f_Pint%0.2f_CI%d_Tmem%d_%ld.dat",L,VIZ,gama,pext,pint,INIT,T_mem,id);
	#endif
	fdensXpcp = fopen(name,"w");


	fprintf(fdensXpcp,"#1:pcp 2:dens00 3:dens01 4:dens10 5:dens11\n");
	dp=0.01;
	pcp=0.;
	pmax=1.;
	while(pcp<=pmax){
		initialState();
		evolve(tmax);
		fprintf(fdensXpcp,"%f ",pcp);
		fflush(fdensXpcp);
		printf("%f ",pcp);
		for(i=0; i<4; ++i){
			fprintf(fdensXpcp,"%f ",dens2[i]);
			fflush(fdensXpcp);
			printf("%f ",dens2[i]);
		}
		fprintf(fdensXpcp,"\n");
		fflush(fdensXpcp);
		printf("\n");
		pcp+=dp;
	}


	fclose(fdensXpcp);
	return;
}
/****************************************************
*   correlation between of 1's in the groups        *
*   as a function of the imitation probability pcp  *
*****************************************************/
void corrXpcp(int tmax){
	int i,j,k,ok,id1,id2,viz,qual1,qual2,smean_e,smean_i,s0e,s0i,sre,sri;
	double dp,pmax;
	unsigned long id;
	char name[100];

	id = (unsigned long)time(NULL);

	ok=0;
	while(ok==0){
	#if  (REDE==0)
		sprintf(name,"corrXpcp_L%d_Gsize%d_G%0.2f_Pext%0.2f_Pint%0.2f_CI%d_Tmem%d_%ld.dat",L,gsize,gama,pext,pint,INIT,T_mem,id);
	#else
		sprintf(name,"corrXpcp_L%d_VIZ%d_G%0.2f_Pext%0.2f_Pint%0.2f_CI%d_Tmem%d_%ld.dat",L,VIZ,gama,pext,pint,INIT,T_mem,id);
	#endif
		fcorrXpcp = fopen(name,"r");
		if(fcorrXpcp!=NULL){
			++id;
			fclose(fcorrXpcp);
		}else{
			ok=1;
		}
	}

	#if  (REDE==0)
	sprintf(name,"corrXpcp_L%d_Gsize%d_G%0.2f_Pext%0.2f_Pint%0.2f_CI%d_Tmem%d_%ld.dat",L,gsize,gama,pext,pint,INIT,T_mem,id);
	#else
	sprintf(name,"corrXpcp_L%d_VIZ%d_G%0.2f_Pext%0.2f_Pint%0.2f_CI%d_Tmem%d_%ld.dat",L,VIZ,gama,pext,pint,INIT,T_mem,id);
	#endif
	fcorrXpcp = fopen(name,"w");


	fprintf(fcorrXpcp,"#1:pcp 3:size_corr 3:corr_i 4:corr_e 5:corr00 6:corr01 7:corr10 8:corr11 9:dens00 10:dens01 11:dens10 12:dens11\n");

	dp=0.01;
	pcp=0.;
	pmax=1.;

	corr_e=(double *)malloc(sizeof(double)*SIZE_CORR+1);
	corr_i=(double *)malloc(sizeof(double)*SIZE_CORR+1);
	for(i=0; i<4; ++i){
		corr_ie[i]=(double *)malloc(sizeof(double)*SIZE_CORR+1);
	}

	while(pcp<=pmax){
		initialState();
		evolve(tmax);

		for(i=0; i<4; ++i){
			memset(corr_ie[i],0.,sizeof(double)*SIZE_CORR+1);//corr_ie[id][distance], id: 1=00, 2=01, 3=10, 4=11
		}
		memset(corr_e,0.,sizeof(double)*SIZE_CORR+1);
		memset(corr_i,0.,sizeof(double)*SIZE_CORR+1);

		corr_e[0]=1.;
		corr_i[0]=1.;
		for(i=0; i<4; ++i){
			corr_ie[i][0]=1.;
		}

		for(i=0; i<N; ++i){
			qual1=i;
			qual2=i;
			s0e=2*s[0][i]-1;
			s0i=2*s[1][i]-1;
			for(j=1; j<=SIZE_CORR; ++j){
				qual1=neighbor[2][qual1];//right
				qual2=neighbor[3][qual2];//down
			//---right neighbor at a distance j---			 
				sre=2*s[0][qual1]-1;
				sri=2*s[1][qual1]-1;
				corr_e[j]+=s0e*sre; 
				corr_i[j]+=s0i*sri; 
			//---bottom neighbor at a distance j---			 
				sre=2*s[0][qual2]-1;
				sri=2*s[1][qual2]-1;
				corr_e[j]+=s0e*sre; 
				corr_i[j]+=s0i*sri; 
			}
		}


		
		for(i=1; i<=SIZE_CORR; ++i){
			corr_e[i]/=(2.*N);//factor 2 is here because corr[i] was counted twice (right and bottom neighbors)
			corr_i[i]/=(2.*N);
		}
	
		smean_e=0.;
		smean_i=0.;
		for(i=0; i<N; ++i){
			smean_e+=2*s[0][i]-1;//turning s=0,1 into s=-1,1
			smean_i+=2*s[1][i]-1;//turning s=0,1 into s=-1,1
		}

		smean_e/=(double)N;
		smean_i/=(double)N;
		
		fprintf(fcorrXpcp,"%f %d %f %f %f %f %f %f\n",pcp,0,corr_i[0],corr_e[0],dens2[0],dens2[1],dens2[2],dens2[3]);
		fflush(fcorrXpcp);
		printf("%f %d %f %f %f %f %f %f\n",pcp,0,corr_i[0],corr_e[0],dens2[0],dens2[1],dens2[2],dens2[3]);
		
		for(i=1; i<=SIZE_CORR;++i){
			corr_e[i]=corr_e[i]-(smean_e*smean_e);
			corr_i[i]=corr_i[i]-(smean_i*smean_i);
			fprintf(fcorrXpcp,"%f %d %f %f ",pcp,i,corr_i[i],corr_e[i]);
			fflush(fcorrXpcp);
			printf("%f %d %f %f ",pcp,i,corr_i[i],corr_e[i]);
			for(j=0; j<4; ++j){
				fprintf(fcorrXpcp,"%f ",dens2[j]);
				fflush(fcorrXpcp);
				printf("%f ",dens2[j]);
			}
			fprintf(fcorrXpcp,"\n");
			printf("\n");
			fflush(fcorrXpcp);
		}
		fprintf(fcorrXpcp,"\n");
		printf("\n");
		fflush(fcorrXpcp);
		
		pcp+=dp;
	}

	//free allocated memory and close open file
	for(i=0; i<4; ++i){
		free(corr_ie[i]);
	}
	free(corr_e);
	free(corr_i);
	fclose(fcorrXpcp);
	return;
}
/************************************************
 *      densities as functions of gamma         * 
 *      factor                                  *
 ************************************************/
void densXgama(int tmax){
	int i,ok;
	double dg,gmax;
	unsigned long id;
	char name[100];

	id = (unsigned long)time(NULL);

	ok=0;
	while(ok==0){
	#if  (REDE==0)
		sprintf(name,"densXgama_L%d_Gsize%d_Pcp%0.2f_Pext%0.2f_Pint%0.2f_Tmem%d_%ld.dat",L,gsize,pcp,pext,pint,T_mem,id);
	#else
		sprintf(name,"densXgama_L%d_VIZ%d_Pcp%0.2f_Pext%0.2f_Pint%0.2f_Tmem%d_%ld.dat",L,VIZ,pcp,pext,pint,T_mem,id);
	#endif
		fdensXgama = fopen(name,"r");
		if(fdensXgama!=NULL){
			++id;
			fclose(fdensXgama);
		}else{
			ok=1;
		}
	}

	#if  (REDE==0)
	sprintf(name,"densXgama_L%d_Gsize%d_Pcp%0.2f_Pext%0.2f_Pint%0.2f_Tmem%d_%ld.dat",L,gsize,pcp,pext,pint,T_mem,id);
	#else
	sprintf(name,"densXgama_L%d_VIZ%d_Pcp%0.2f_Pext%0.2f_Pint%0.2f_Tmem%d_%ld.dat",L,VIZ,pcp,pext,pint,T_mem,id);
	#endif
	fdensXgama = fopen(name,"w");


	fprintf(fdensXgama,"#1:gama 2:dens00 3:dens01 4:dens10 5:dens11\n");
	dg=0.01;
	gama=0.;
	gmax=1.;
	while(gama<=gmax){
		initialState();
		evolve(tmax);
		fprintf(fdensXgama,"%f ",gama);
		fflush(fdensXgama);
		printf("%f ",gama);
		for(i=0; i<4; ++i){
			fprintf(fdensXgama,"%f ",(double)dens2[i]);
			fflush(fdensXgama);
			printf("%f ",(double)dens2[i]);
		}
		fprintf(fdensXgama,"\n");
		fflush(fdensXgama);
		printf("\n");
		gama+=dg;
	}


	fclose(fdensXgama);
	return;
}
/************************************************************
 *      diagram for the densities as functions of the       *
 *      imitation and internalization probabilities         *
 ************************************************************/
void diagram_pcpXpintXd(int tmax){
	int ok;
	double dp, pmax;
	unsigned long id;
	char name[100];

	id = (unsigned long)time(NULL);

	ok=0;
	while(ok==0){
	#if  (REDE==0)
		sprintf(name,"diagram_pcpXpintXd_L%d_Gsize%d_G%0.2f_Pext%0.2f_Tmem%d_%ld.dat",L,gsize,gama,pext,T_mem,id);
	#else
		sprintf(name,"diagram_pcpXpintXd_L%d_VIZ%d_G%0.2f_Pext%0.2f_Tmem%d_%ld.dat",L,VIZ,gama,pext,T_mem,id);
	#endif
		fdiagpcpi = fopen(name,"r");
		if(fdiagpcpi!=NULL){
			++id;
			fclose(fdiagpcpi);
		}else{
			ok=1;
		}
	}

#if  (REDE==0)
	sprintf(name,"diagram_pcpXpintXd_L%d_Gsize%d_G%0.2f_Pext%0.2f_Tmem%d_%ld.dat",L,gsize,gama,pext,T_mem,id);
	#else
	sprintf(name,"diagram_pcpXpintXd_L%d_VIZ%d_G%0.2f_Pext%0.2f_Tmem%d_%ld.dat",L,VIZ,gama,pext,T_mem,id);
	#endif
	fdiagpcpi = fopen(name,"w");

	fprintf(fdiagpcpi,"#1:pcp 2:pint 3:dens00 4:dens01 5:dens10 6:dens11\n");
	dp=0.01;
	pcp=0.;
	pmax=1.;
	while(pcp<=pmax){
		pint=dp;
		while(pint<pmax){
			initialState();
			evolve(tmax);
			fprintf(fdiagpcpi,"%f %f %f %f %f %f\n",pcp,pint,dens2[0],dens2[1],dens2[2],dens2[3]);
			fflush(fdiagpcpi);
			printf("%f %f %f %f %f %f\n",pcp,pint,dens2[0],dens2[1],dens2[2],dens2[3]);
			pint+=dp;
		}
		pcp+=dp;
	}


	fclose(fdiagpcpi);
	return;
}
/************************************************************
 *      diagram for the densities as functions of the       *
 *      imitation prob and gama factor                      *
 ************************************************************/
void diagram_pcpXgamaXd(int tmax){
		int ok;
		double dp, pmax;
		unsigned long id;
		char name[100];

		id = (unsigned long)time(NULL);

		ok=0;
		while(ok==0){
		#if  (REDE==0)
			sprintf(name,"diagram_pcpXgamaXd_L%d_Gsize%d_Pint%0.2f_Pext%0.2f_Tmem%d_%ld.dat",L,gsize,pint,pext,T_mem,id);
		#else
			sprintf(name,"diagram_pcpXgamaXd_L%d_VIZ%d_Pint%0.2f_Pext%0.2f_Tmem%d_%ld.dat",L,VIZ,pint,pext,T_mem,id);
		#endif
			fdiagpgd = fopen(name,"r");
			if(fdiagpgd!=NULL){
				++id;
				fclose(fdiagpgd);
			}else{
				ok=1;
			}
		}

		#if  (REDE==0)
		sprintf(name,"diagram_pcpXgamaXd_L%d_Gsize%d_Pint%0.2f_Pext%0.2f_Tmem%d_%ld.dat",L,gsize,pint,pext,T_mem,id);
		#else
		sprintf(name,"diagram_pcpXgamaXd_L%d_VIZ%d_Pint%0.2f_Pext%0.2f_Tmem%d_%ld.dat",L,VIZ,pint,pext,T_mem,id);
		#endif
		fdiagpgd = fopen(name,"w");


		fprintf(fdiagpgd,"#1:pcp 2:gama  3:dens00 4:dens01 5:dens10 6:dens11\n");
		dp=0.01;
		pcp=0.;
		pmax=1.;
		while(pcp<=pmax){
			gama=dp;
			while(gama<pmax){
				initialState();
				evolve(tmax);
				fprintf(fdiagpgd,"%f %f %f %f %f %f\n",pcp,gama,dens2[0],dens2[1],dens2[2],dens2[3]);
				fflush(fdiagpgd);
				printf("%f %f %f %f %f %f\n",pcp,gama,dens2[0],dens2[1],dens2[2],dens2[3]);
				gama+=dp;
			}
			pcp+=dp;
		}


		fclose(fdiagpgd);
		return;
	}
/***************************************************************
*        Stores spatial configurations and creates             *
*        scripts for gnuplot                                   *
****************************************************************/
void save_config(int t){
	int i,j,g,ok,idj,k,n;
	double pointsize=3.5,*averoverlap;
	char name[150],namedat[180],*name_gpe,*name_gpi,*name_gphd,*name_gp4s;
	FILE *fconfig,*fgp_ext,*fgp_int,*fgp_hXd,*fgp_4s;

	name_gpe=(char *)malloc(sizeof(char)*180);//name of the files for external states
	name_gpi=(char *)malloc(sizeof(char)*180);//name of the files for internal states
	name_gphd=(char *)malloc(sizeof(char)*180);//name of the files for harmonial versus desarmonious states
	name_gp4s=(char *)malloc(sizeof(char)*180);//name of the files for the 4 combination of states (00, 01, 10 and 11)

#if  (REDE==0)
	sprintf(name,"snapshot_L%d_Gsize%d_G%0.2f_Pcp%0.2f_Pext%0.2f_Pint%0.2f_CI%d_Di1%0.2f_Tmem%d",L,gsize,gama,pcp,pext,pint,INIT,Di1,T_mem);
#else
	sprintf(name,"snapshot_L%d_VIZ%d_G%0.2f_Pcp%0.2f_Pext%0.2f_Pint%0.2f_CI%d_Di1%0.2f_Tmem%d",L,VIZ,gama,pcp,pext,pint,INIT,Di1,T_mem);
#endif

#if (SAVE_CONFIG_ID==0)//files for making animations (they have to be in an ascending order)
	sprintf(namedat,"%s_idt%f.dat",name,(double)t/Teq);
	sprintf(name_gpe,"%s-ext_idt%f.gp",name,(double)t/Teq);
	sprintf(name_gpi,"%s-int_idt%f.gp",name,(double)t/Teq);
	sprintf(name_gphd,"%s-harmXdesarm_idt%f.gp",name,(double)t/Teq);
	sprintf(name_gp4s,"%s_idt%f.gp",name,(double)t/Teq);
#endif
#if (SAVE_CONFIG_ID==1)
	sprintf(namedat,"%s_Tf%d.dat",name,t);
	sprintf(name_gpe,"%s-ext_Tf%d.gp",name,t);
	sprintf(name_gpi,"%s-int_Tf%d.gp",name,t);
	sprintf(name_gphd,"%s-harmXdesarm_Tf%d.gp",name,t);
	sprintf(name_gp4s,"%s_Tf%d.gp",name,t);
#endif

	fconfig = fopen(namedat,"w");
	fgp_ext = fopen(name_gpe,"w");
	fgp_int = fopen(name_gpi,"w");
	fgp_hXd = fopen(name_gphd,"w");
	fgp_4s = fopen(name_gp4s,"w");


	//storing states for each x and y positions (state 0=00, state 1=01, state 2=10, stae 3=11)
	for(i = 0; i < L; ++i){//coluna
		for(j = 0; j < L; ++j){//linha
			//#1:x 2:y 3:state(from 0 to 4) 4:internal-state  5:external-state
			fprintf(fconfig,"%d %d %d %d %d\n",i,j,2*s[0][i+j*L]+s[1][i+j*L],s[0][i+j*L],s[1][i+j*L]);
		}
		fprintf(fconfig,"\n");
	}
/*-------------creating scripts for gnuplot images------------------*/
#if (SAVE_CONFIG_ID==0)//png, for creating animations
	fprintf(fgp_ext,"set term png size 720,540\n");
	fprintf(fgp_ext,"set output'%s-ext_idt%f.png'\n",name,(double)t/Teq);

	fprintf(fgp_int,"set term png size 720,540\n");
	fprintf(fgp_int,"set output'%s-int_idt%f.png'\n",name,(double)t/Teq);

	fprintf(fgp_hXd,"set term png size 720,540\n");
	fprintf(fgp_hXd,"set output'%s-harmXdesarm_idt%f.png'\n",name,(double)t/Teq);

	fprintf(fgp_4s,"set term png size 720,540\n");
	fprintf(fgp_4s,"set output'%s_idt%f.png'\n",name,(double)t/Teq);
#endif
#if (SAVE_CONFIG_ID==1)//eps for creating single images
	fprintf(fgp_ext,"set term post eps enha color 20\n");
	fprintf(fgp_ext,"set output'%s-ext_Tf%d_fm.eps'\n",name,t);

	fprintf(fgp_int,"set term post eps enha color 20\n");
	fprintf(fgp_int,"set output'%s-int_Tf%d_fm.eps'\n",name,t);

	fprintf(fgp_hXd,"set term post eps enha color 20\n");
	fprintf(fgp_hXd,"set output'%s-harmXdesarm_Tf%d_fm.eps'\n",name,t);

	fprintf(fgp_4s,"set term post eps enha color 20\n");
	fprintf(fgp_4s,"set output'%s_Tf%d_fm.eps'\n",name,t);
#endif

	/*external state files*/
	fprintf(fgp_ext,"\n set title'{/=15 di=%0.2f, de=%0.2f, di(0)=%0.2f, CI=%d}'\n",dens[0],dens[1],Di1,INIT);
	fprintf(fgp_ext,"\n set key outside\n");
	fprintf(fgp_ext,"\n set key title'{/=12 External}'\n");
	fprintf(fgp_ext,"\n set key box\n");
	fprintf(fgp_ext,"unset xtics\n");
	fprintf(fgp_ext,"unset ytics\n");
	fprintf(fgp_ext,"L=%d\n",L);
	fprintf(fgp_ext,"set yr[0:L-1]\n");
	fprintf(fgp_ext,"set xr[0:L-1]\n");
	fprintf(fgp_ext,"set size square\n");
	fprintf(fgp_ext,"set pointsize %0.2lf\n",pointsize);
	fprintf(fgp_ext,"plot '%s' u 1:(($5==0)?$2:1/0) w p pt 5 lc rgb 'red' t'{/=12 0}',\\\n",namedat);
	fprintf(fgp_ext,"'%s' u 1:(($5==1)?$2:1/0) w p pt 5 lc rgb 'blue' t'{/=12 1}',\\\n",namedat);

	/*internal state files*/
	fprintf(fgp_int,"\n set title'{/=15 di=%0.2f, de=%0.2f, di(0)=%0.2f, CI=%d}'\n",dens[0],dens[1],Di1,INIT);
	fprintf(fgp_int,"\n set key outside\n");
	fprintf(fgp_int,"\n set key title'{/=12 Internal}'\n");
	fprintf(fgp_int,"\n set key box\n");
	fprintf(fgp_int,"unset xtics\n");
	fprintf(fgp_int,"unset ytics\n");
	fprintf(fgp_int,"L=%d\n",L);
	fprintf(fgp_int,"set yr[0:L-1]\n");
	fprintf(fgp_int,"set xr[0:L-1]\n");
	fprintf(fgp_int,"set size square\n");
	fprintf(fgp_int,"set pointsize %0.2lf\n",pointsize);
	fprintf(fgp_int,"plot '%s' u 1:(($4==0)?$2:1/0) w p pt 5 lc rgb 'red' t'{/=12 0}',\\\n",namedat);
	fprintf(fgp_int,"'%s' u 1:(($4==1)?$2:1/0) w p pt 5 lc rgb 'blue' t'{/=12 1}',\\\n",namedat);
	
	/*harmonial versus desarmonial states files*/
	fprintf(fgp_hXd,"\n set title'{/=15 di=%0.2f, de=%0.2f,di(0)=%0.2f, CI=%d}'\n",dens[0],dens[1], Di1,INIT);
	fprintf(fgp_hXd,"\n set key outside\n");
	fprintf(fgp_hXd,"\n set key title'{/=12 Harm. Vs Desarm.}'\n");
	fprintf(fgp_hXd,"\n set key box\n");
	fprintf(fgp_hXd,"unset xtics\n");
	fprintf(fgp_hXd,"unset ytics\n");
	fprintf(fgp_hXd,"L=%d\n",L);
	fprintf(fgp_hXd,"set yr[0:L-1]\n");
	fprintf(fgp_hXd,"set xr[0:L-1]\n");
	fprintf(fgp_hXd,"set size square\n");
	fprintf(fgp_hXd,"set pointsize %0.2lf\n",pointsize);
	fprintf(fgp_hXd,"plot '%s' u 1:((($3==1)||($3==2))?$2:1/0) w p pt 5 lc rgb 'red' t'{/=12 01+10}',\\\n",namedat);
	fprintf(fgp_hXd,"'%s' u 1:((($3==0)||($3==3))?$2:1/0) w p pt 5 lc rgb 'blue' t'{/=12 00+11}',\\\n",namedat);

	/*4 combinations of states files*/
	fprintf(fgp_4s,"\n set title'{/=15 di=%0.2f, de=%0.2f,di(0)=%0.2f, CI=%d}'\n",dens[0],dens[1], Di1,INIT);
	fprintf(fgp_4s,"\n set key outside\n");
	fprintf(fgp_4s,"\n set key title'{/=12 4 states}'\n");
	fprintf(fgp_4s,"\n set key box\n");
	fprintf(fgp_4s,"unset xtics\n");
	fprintf(fgp_4s,"unset ytics\n");
	fprintf(fgp_4s,"L=%d\n",L);
	fprintf(fgp_4s,"set yr[0:L-1]\n");
	fprintf(fgp_4s,"set xr[0:L-1]\n");
	fprintf(fgp_4s,"set size square\n");
	fprintf(fgp_4s,"set pointsize %0.2lf\n",pointsize);
	fprintf(fgp_4s,"plot '%s' u 1:(($3==0)?$2:1/0) w p pt 5 lc rgb 'red' t'{/=12 00}',\\\n",namedat);
	fprintf(fgp_4s,"'%s' u 1:(($3==1)?$2:1/0) w p pt 5 lc rgb 'green' t'{/=12 01}',\\\n",namedat);
	fprintf(fgp_4s,"'%s' u 1:(($3==2)?$2:1/0) w p pt 5 lc rgb '#FF00FF' t'{/=12 10}',\\\n",namedat);
	fprintf(fgp_4s,"'%s' u 1:(($3==3)?$2:1/0) w p pt 5 lc rgb 'blue' t'{/=12 11}'\n",namedat);

/*-----free allocated memory and close open files0-----------*/

       free(averoverlap);
       free(name_gpe);
       free(name_gpi);
       free(name_gphd);
       free(name_gp4s);
       fclose(fconfig);
       fclose(fgp_ext);
       fclose(fgp_int);
       fclose(fgp_hXd);
       fclose(fgp_4s);

/*---------------------------*/
       return;
}
/**************************************
 *  free allocated memory and         *
 *  close files that are still open   *
 *************************************/
void freeMemory(void){

	free(s[0]);
        free(s[1]);
	free(count_mem);
	free(mem);
	free(mem0);
#if  (REDE==1)
        free(neighbor);
#endif
#if  (REDE==0)
	free(list);
	free(group);
#endif

        return;
}


