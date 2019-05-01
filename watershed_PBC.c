// gcc -O3 -o watershed watershed.c
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <stdbool.h>

struct NEIGHBORHOOD
{
    int64_t NNZV; // Compressed Row "Binary" Storage (CRBS)
    int64_t *IV; // Compressed Row "Binary" Storage (CRBS)
    int32_t *JV; // Compressed Row "Binary" Storage (CRBS)
};

struct LATTICE
{
    int32_t lx; // Lattice size on X dimension
    int32_t ly; // Lattice size on Y dimension
    int32_t n; // Total number of sites
    double lon; // Longitude of the most northwest site
    double lat; // Latitude of the most northwest site
    double delta; // increment in longitude/latitude
    double R; // Celestial body's radius
    
    float *h; // Heights
    
    struct NEIGHBORHOOD *nei; // Lattice neighboorhood
};

struct LIST
{
    int32_t n; // Length
    int32_t m; // Size
    int32_t *key; // Keys
};

struct HEAP
{
    int32_t n; // Length
    int32_t m; // Size
    float *key; // Keys
    int32_t *item; // Positions
};

struct SYSTEM
{
    struct LATTICE *L; // Lattice
    
    int32_t *sigma; // Basins
    struct LIST *sink; // Sinks
    
    struct LIST *stalker; // Stalkers
    int32_t *check; // Check
    
    struct LIST *burner; // Burners
    int32_t *status; // Status
    
    struct HEAP *H; // Binary Heap
};

void input(char name[100],struct SYSTEM *S);
void set_neighbors(struct SYSTEM *S);
void coastline(char hcutoff[100],struct SYSTEM *S);
void drainage_basin(struct SYSTEM *S);
void invasion_percolation(struct SYSTEM *S,int32_t k);
void output(char name[100],char hcutoff[100],struct SYSTEM *S);

int32_t parent(int32_t i);
int32_t left(int32_t i);
int32_t right(int32_t i);
void min_heap_insert(struct HEAP *A,float key,int32_t item);
void heap_decrease_key(struct HEAP *A,int32_t i,float key,int32_t item);
void heap_extract_min(struct HEAP *A,float *key,int32_t *item);
void min_heapify(struct HEAP *A,int32_t i);

void set_free(struct SYSTEM *S);

int32_t main(int32_t argc,char *argv[])
{
    char name[100];
    char hcutoff[100];
    
    struct SYSTEM *S;
    
    strcpy(name,argv[1]);
    strcpy(hcutoff,argv[2]);
    
    S=(struct SYSTEM *)malloc(sizeof(struct SYSTEM));
    input(name,S);
    coastline(hcutoff,S);
    drainage_basin(S);
    output(name,hcutoff,S);
    set_free(S);
    
    return(0);
}

void input(char name[100],struct SYSTEM *S)
{
    char filename[400];
    FILE *f;
    
    int32_t lx;
    int32_t ly;
    double lon;
    double lat;
    double delta;
    double R;
    
    int32_t i,j,k;
    float height;
    
    int32_t t;
    
    strcpy(filename,name);
    strcat(filename,".dat");
    f=fopen(filename,"r");
    fscanf(f,"%d %d %lf %lf %lf %lf\n",&lx,&ly,&lon,&lat,&delta,&R);
    
    S->L=(struct LATTICE *)malloc(sizeof(struct LATTICE));
    S->L->lx=lx;
    S->L->ly=ly;
    S->L->n=S->L->lx*S->L->ly;
    S->L->lon=lon;
    S->L->lat=lat;
    S->L->delta=delta;
    S->L->R=R;
    S->L->h=(float *)malloc(S->L->n*sizeof(float));
    
    set_neighbors(S); // Neighborhood
    
    srand(123); // Seed
    for(t=0;t<S->L->n;t++)
    {
        fscanf(f,"%d %d %f\n",&i,&j,&height);
        
        k=i+j*S->L->lx;
        
        height=height+0.001*(float)rand()/(float)RAND_MAX; // Noise
        
        S->L->h[k]=height;
    }
    fclose(f);
}

void set_neighbors(struct SYSTEM *S)
{
    int32_t i,j,k;
    int32_t r,s;
    int64_t t;
    
    /*PBC*/
    S->L->nei=(struct NEIGHBORHOOD *)malloc(sizeof(struct NEIGHBORHOOD));
    
    S->L->nei->NNZV=2*(int64_t)S->L->lx*((int64_t)S->L->lx+2*(int64_t)S->L->ly-4);
    S->L->nei->IV=(int64_t *)malloc((S->L->n+1)*sizeof(int64_t));
    S->L->nei->JV=(int32_t *)malloc(S->L->nei->NNZV*sizeof(int32_t));
    /*PBC*/
    
    /*PBC*/
    t=0;
    S->L->nei->IV[0]=t;
    for(k=0;k<S->L->n;k++)
    {
        i=k%S->L->lx;
        j=k/S->L->lx;
        
        if(j==0)
        {
            if(i==0)
            {
                /*VON NEUMMAN*/
                S->L->nei->JV[t]=k+1;
                t=t+1;
                for(r=0;r<S->L->lx;r++)
				{
					s=r;
					
					if(s!=k+S->L->lx-1 && s!=k && s!=k+1)
					{
						S->L->nei->JV[t]=s; // PBC
						t=t+1;
					}
				}
				S->L->nei->JV[t]=k+S->L->lx-1; // PBC
				t=t+1;
                S->L->nei->JV[t]=k+S->L->lx;
                t=t+1;
                
                S->L->nei->IV[k+1]=t;
                /*VON NEUMMAN*/
            }
            else
            {
                if(i==S->L->lx-1)
                {
                    /*VON NEUMMAN*/
                    S->L->nei->JV[t]=k-S->L->lx+1; // PBC
					t=t+1;
					for(r=0;r<S->L->lx;r++)
					{
						s=r;
						
						if(s!=k-1 && s!=k && s!=k-S->L->lx+1)
						{
							S->L->nei->JV[t]=s; // PBC
							t=t+1;
						}
					}
                    S->L->nei->JV[t]=k-1;
                    t=t+1;
                    S->L->nei->JV[t]=k+S->L->lx;
                    t=t+1;
                    
                    S->L->nei->IV[k+1]=t;
                    /*VON NEUMMAN*/
                }
                else
                {
                    /*VON NEUMMAN*/
                    for(r=0;r<S->L->lx;r++)
					{
						s=r;
						
						if(s!=k-1 && s!=k && s!=k+1)
						{
							S->L->nei->JV[t]=s; // PBC
							t=t+1;
						}
						else
						{
							if(s==k-1)
							{
								S->L->nei->JV[t]=s;
								t=t+1;
							}
							
							if(s==k+1)
							{
								S->L->nei->JV[t]=s;
								t=t+1;
							}
						}
					}
                    S->L->nei->JV[t]=k+S->L->lx;
                    t=t+1;
                    
                    S->L->nei->IV[k+1]=t;
                    /*VON NEUMMAN*/
                }
            }
        }
        else
        {
            if(j<S->L->ly-1) //k<=(S->L->ly-1)*S->L->lx-1
            {
                if(i==0)
                {
                    /*VON NEUMMAN*/
                    S->L->nei->JV[t]=k-S->L->lx;
                    t=t+1;
                    S->L->nei->JV[t]=k+1;
                    t=t+1;
                    S->L->nei->JV[t]=k-1+S->L->lx; // PBC
					t=t+1;
                    S->L->nei->JV[t]=k+S->L->lx;
                    t=t+1;
                    
                    S->L->nei->IV[k+1]=t;
                    /*VON NEUMMAN*/
                }
                else
                {
                    if(i==S->L->lx-1)
                    {
                        /*VON NEUMMAN*/
                        S->L->nei->JV[t]=k-S->L->lx;
                        t=t+1;
                        S->L->nei->JV[t]=k-S->L->lx+1; // PBC
						t=t+1;
                        S->L->nei->JV[t]=k-1;
                        t=t+1;
                        S->L->nei->JV[t]=k+S->L->lx;
                        t=t+1;
                        
                        S->L->nei->IV[k+1]=t;
                        /*VON NEUMMAN*/
                    }
                    else
                    {
                        /*VON NEUMMAN*/
                        S->L->nei->JV[t]=k-S->L->lx;
                        t=t+1;
                        S->L->nei->JV[t]=k-1;
                        t=t+1;
                        S->L->nei->JV[t]=k+1;
                        t=t+1;
                        S->L->nei->JV[t]=k+S->L->lx;
                        t=t+1;
                        
                        S->L->nei->IV[k+1]=t;
                        /*VON NEUMMAN*/
                    }
                }
            }
            else
            {
                if(i==0)
                {
                    /*VON NEUMMAN*/
                    S->L->nei->JV[t]=k-S->L->lx;
                    t=t+1;
                    S->L->nei->JV[t]=k+1;
                    t=t+1;
                    for(r=0;r<S->L->lx;r++)
					{
						s=r+(S->L->ly-1)*S->L->lx;
						
						if(s!=k-1+S->L->lx && s!=k && s!=k+1)
						{
							S->L->nei->JV[t]=s; // PBC
							t=t+1;
						}
					}
					S->L->nei->JV[t]=k-1+S->L->lx; // PBC
					t=t+1;
                    
                    S->L->nei->IV[k+1]=t;
                    /*VON NEUMMAN*/
                }
                else
                {
                    if(i==S->L->lx-1)
                    {
                        /*VON NEUMMAN*/
                        S->L->nei->JV[t]=k-S->L->lx;
                        t=t+1;
                        S->L->nei->JV[t]=k-S->L->lx+1; // PBC
						t=t+1;						
						for(r=0;r<S->L->lx;r++)
						{
							s=r+(S->L->ly-1)*S->L->lx;
							
							if(s!=k-1 && s!=k && s!=k-S->L->lx+1)
							{
								S->L->nei->JV[t]=s; // PBC
								t=t+1;
							}
						}
                        S->L->nei->JV[t]=k-1;
                        t=t+1;
                        
                        S->L->nei->IV[k+1]=t;
                        /*VON NEUMMAN*/
                    }
                    else
                    {
                        /*VON NEUMMAN*/
                        S->L->nei->JV[t]=k-S->L->lx;
                        t=t+1;
                        for(r=0;r<S->L->lx;r++)
						{
							s=r+(S->L->ly-1)*S->L->lx;
							
							if(s!=k-1 && s!=k && s!=k+1)
							{
								S->L->nei->JV[t]=s; // PBC
								t=t+1;
							}
							else
							{
								if(s==k-1)
								{
									S->L->nei->JV[t]=s;
									t=t+1;
								}
								
								if(s==k+1)
								{
									S->L->nei->JV[t]=s;
									t=t+1;
								}
							}
						}
                        
                        S->L->nei->IV[k+1]=t;
                        /*VON NEUMMAN*/
                    }
                }
            }
        }
    }
    /*PBC*/
}

void coastline(char hcutoff[100],struct SYSTEM *S)
{
    int32_t k;
    int64_t r;
    int32_t s;
    
    S->sigma=(int32_t *)malloc(S->L->n*sizeof(int32_t));
    
    S->sink=(struct LIST *)malloc(sizeof(struct LIST));
    S->sink->n=S->L->n;
    S->sink->m=0;
    S->sink->key=(int32_t *)malloc(S->sink->n*sizeof(int32_t));
    
    for(k=0;k<S->L->n;k++)
    {
        if(S->L->h[k]>atof(hcutoff))
		{
		    S->sigma[k]=-1;
		    
            for(r=S->L->nei->IV[k];r<S->L->nei->IV[k+1];r++) // from IA[i] to IA[i+1]-1
            {
                s=S->L->nei->JV[r];
                
                if(S->L->h[s]<=atof(hcutoff))
                {
                    S->sink->key[S->sink->m]=k;
                    S->sink->m++;
                    
                    break;
                }
            }
        }
        else
        {
            S->sigma[k]=-2;
        }
    }
}

void drainage_basin(struct SYSTEM *S)
{
    int32_t q;
    
    int32_t i,j,k;
    int64_t r,u;
    int32_t s,t;
    int32_t v;
    
    S->stalker=(struct LIST *)malloc(sizeof(struct LIST));
    S->stalker->n=S->L->n;
    S->stalker->m=0;
    S->stalker->key=(int32_t *)malloc(S->stalker->n*sizeof(int32_t));
    
    S->check=(int32_t *)malloc(S->L->n*sizeof(int32_t));
    
    S->burner=(struct LIST *)malloc(sizeof(struct LIST));
    S->burner->n=S->L->n;
    S->burner->m=0;
    S->burner->key=(int32_t *)malloc(S->burner->n*sizeof(int32_t));
    
    S->status=(int32_t *)malloc(S->L->n*sizeof(int32_t));
    
    S->H=(struct HEAP *)malloc(sizeof(struct HEAP));
    S->H->n=S->L->n;
    //S->H->m=0;
    S->H->key=(float *)malloc(S->H->n*sizeof(float));
    S->H->item=(int32_t *)malloc(S->H->n*sizeof(int32_t));
    
    for(i=0;i<S->L->n;i++)
    {
        S->check[i]=-2;
    }
    
    for(i=0;i<S->L->n;i++)
    {
        if(S->sigma[i]==-2)
        {
            S->status[i]=-2;
        }
        else
        {
            S->status[i]=-1;
        }
    }
    
    for(i=0;i<S->sink->m;i++)
    {
        k=S->sink->key[i];
        
        invasion_percolation(S,k); // S->sigma[k]=k
        
        for(j=0;j<S->stalker->m;j++)
        {
            t=S->stalker->key[j];
            S->check[t]=-2;
        }
        
        S->stalker->m=0; // Restarting
        
        S->stalker->key[S->stalker->m]=k;
        S->stalker->m++;
        S->check[k]=1;
        
        q=0;
        while(q<S->stalker->m)
        {
            t=S->stalker->key[q];
            
            for(r=S->L->nei->IV[t];r<S->L->nei->IV[t+1];r++) // from IA[i] to IA[i+1]-1
            {
                s=S->L->nei->JV[r];
                
                if(S->sigma[s]>=-1)
                {
                    invasion_percolation(S,s);
                    
                    for(u=S->L->nei->IV[s];u<S->L->nei->IV[s+1];u++) // from IA[i] to IA[i+1]-1
                    {
                        v=S->L->nei->JV[u];
                        
                        if(S->sigma[v]>=-1)
                        {
                            invasion_percolation(S,v);
                            
                            if(S->sigma[s]==S->sigma[k])
                            {
                                if(S->sigma[v]!=S->sigma[k])
                                {
                                    if(S->check[s]==-2)
                                    {
                                        S->stalker->key[S->stalker->m]=s;
                                        S->stalker->m++;
                                        S->check[s]=1;
                                    }
                                }
                            }
                            else
                            {
                                if(S->sigma[v]==S->sigma[k])
                                {
                                    if(S->check[v]==-2)
                                    {
                                        S->stalker->key[S->stalker->m]=v;
                                        S->stalker->m++;
                                        S->check[v]=1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            q++;
        }
    }
}

void invasion_percolation(struct SYSTEM *S,int32_t k)
{
    bool stop;
    
    float key;
    int32_t item;
    
    int32_t i;
    int64_t r;
    int32_t s,t;
    int32_t w;
    
    if(S->sigma[k]==-1)
    {
        for(t=0;t<S->burner->m;t++)
        {
            i=S->burner->key[t];
            S->status[i]=-1;
        }
        
        S->burner->m=0; // Restarting
        S->H->m=0; // Restarting
        
        S->burner->key[S->burner->m]=k;
        S->burner->m++;
        S->status[k]=1;
        
        min_heap_insert(S->H,S->L->h[k],k);
        
        stop=false;
        while(S->H->m>0 && stop==false)
        {
            heap_extract_min(S->H,&key,&item);
            
            t=item;
            
            S->status[t]=2;
            
            for(r=S->L->nei->IV[t];r<S->L->nei->IV[t+1];r++) // from IA[i] to IA[i+1]-1
            {
                s=S->L->nei->JV[r];
                
                if(S->status[s]==-2)
                {
                    w=t;
                    stop=true;
                    
                    break;
                }
            }
            
            if(stop==false)
            {
                for(r=S->L->nei->IV[t];r<S->L->nei->IV[t+1];r++) // from IA[i] to IA[i+1]-1
                {
                    s=S->L->nei->JV[r];
                    
                    if(S->status[s]==-1)
                    {
                        S->burner->key[S->burner->m]=s;
                        S->burner->m++;
                        S->status[s]=1;
                        min_heap_insert(S->H,S->L->h[s],s);
                    }
                }
            }
        }
        
        if(stop==true)
        {
            S->sigma[k]=w;
        }
    }
}

void output(char name[100],char hcutoff[100],struct SYSTEM *S)
{
    FILE *f;
    char filename[400];
    
    int32_t i,j,k;
    
    strcpy(filename,"sigma_");
    strcat(filename,name);
    strcat(filename,"_h"); // hcutoff
    strcat(filename,hcutoff); // hcutoff
    strcat(filename,".dat");
    f=fopen(filename,"w");
    fprintf(f,"%d %d %lf %lf %lf %lf\n",S->L->lx,S->L->ly,S->L->lon,S->L->lat,S->L->delta,S->L->R);
    for(k=0;k<S->L->n;k++)
    {
        i=k%S->L->lx;
        j=k/S->L->lx;
        
        fprintf(f,"%d %d %d\n",i,j,S->sigma[k]);
    }
    fclose(f);
}

int32_t parent(int32_t i)
{
    return (i-1)/2;
}

int32_t left(int32_t i)
{
    return 2*i+1;
}

int32_t right(int32_t i)
{
    return 2*i+2;
}

void min_heap_insert(struct HEAP *A,float key,int32_t item)
{
    A->m++;
    A->key[A->m-1]=FLT_MAX;
    
    heap_decrease_key(A,A->m-1,key,item);
}

void heap_decrease_key(struct HEAP *A, int32_t i, float key,int32_t item)
{
    int32_t p;
    float aux_key;
    int32_t aux_item;
    
    if(key>A->key[i])
    {
        printf("New key is larger than current key\n");
        exit(0);
    }
    
    A->key[i]=key;
    A->item[i]=item;
    
    p=parent(i);
    while(i>0 && A->key[p]>A->key[i])
    {
        aux_key=A->key[i];
        A->key[i]=A->key[p];
        A->key[p]=aux_key;
        
        aux_item=A->item[i];
        A->item[i]=A->item[p];
        A->item[p]=aux_item;
        
        i=p;
        p=parent(i);
    }
}

void heap_extract_min(struct HEAP *A, float *key, int32_t *item)
{
    if(A->m<1)
    {
        printf("Heap underflow\n");
        exit(0);
    }
    
    *key=A->key[0];
    *item=A->item[0];
    
    A->key[0]=A->key[A->m-1];
    A->item[0]=A->item[A->m-1];
    
    A->m--;
    
    min_heapify(A,0);
}

void min_heapify(struct HEAP *A,int32_t i)
{
    int32_t l;
    int32_t r;
    int32_t smallest;
    float aux_key;
    int32_t aux_item;
    
    l=left(i);
    r=right(i);
    
    if(l<A->m && A->key[l]<A->key[i])
    {
        smallest=l;
    }
    else
    {
        smallest=i;
    }
    
    if(r<A->m && A->key[r]<A->key[smallest])
    {
        smallest=r;
    }
    
    if(smallest!=i)
    {
        aux_key=A->key[i];
        A->key[i]=A->key[smallest];
        A->key[smallest]=aux_key;
        
        aux_item=A->item[i];
        A->item[i]=A->item[smallest];
        A->item[smallest]=aux_item;
        
        min_heapify(A,smallest);
    }
}

void set_free(struct SYSTEM *S)
{
    free(S->H->item);
    free(S->H->key);
    free(S->H);
    
    free(S->status);
    free(S->burner->key);
    free(S->burner);
    
    free(S->check);
    free(S->stalker->key);
    free(S->stalker);
    
    free(S->sink->key);
    free(S->sink);
    
    free(S->sigma);
    
    free(S->L->nei->JV);
    free(S->L->nei->IV);
    free(S->L->nei);
    
    free(S->L->h);
    free(S->L);
    free(S);
}
