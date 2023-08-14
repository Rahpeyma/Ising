#include<stdio.h>
#include<time.h>
#include<stdlib.h>
#include<math.h>
#define ROWS 16
#define COLS 16
#define N2 COLS*ROWS
#define five 5
#define NN 900 

int energy(int i,int j,int s[ROWS][COLS]){
   int im1=i-1,jm1=j-1,ip1=i+1,jp1=j+1;
// return -s[i][j]*(s[(i-1+ROWS)%ROWS][j]+s[(i+1)%ROWS][j]+s[i][(j+1)%COLS]+s[i][(j-1+COLS)%COLS]);
    if (i==0) im1=ROWS-1;
    else if (i==ROWS-1) ip1=0;
    if (j==0) jm1=COLS-1;
    else if (j==ROWS-1) jp1=0;
    return -s[i][j]*(s[(im1)][j]+s[(ip1)][j]+s[i][(jp1)]+s[i][(jm1)]);
}
float prob(int energy,float array[five]){
    // return exp(-energy/T);
    float ans;
    switch (energy)
    {
    case 8:
        ans= array[0];
        break;
    case -8:
        ans= array[1];
        break;
    case 4:
        ans= array[2];
        break;
    case -4:
        ans= array[3];
        break;
    case 0:
        ans= array[4];
        break;
    }
    return ans;
}
float M(int s[ROWS][COLS]){
    int sum=0;
    for (int i = 0; i < ROWS; i++)
    {
        for (int j = 0; j < COLS; j++)
        {
            sum+=s[i][j];
        }
    }
    return sum/((float)N2);
}
float energy_t(int s[ROWS][COLS]){
    int sum=0;
    for (int i = 0; i < ROWS; i++)
    {
        for (int j = 0; j < COLS; j++)
        {
            sum+=energy(i,j,s);
        }
    }
    return sum/((float)N2);
}

int main(){
    srand(time(NULL));
    int s[ROWS][COLS],ii;
    float temp[401];
    int x,y,e1;
    float M_main=0.0,MS,arrayexp[five];
    float M2_main=0.0,yrand;
    float energy_main=0.0,energy2_main=0.0,EnergyTotal=0.0;
    float X,Cv;

    for(int i=0;i<ROWS;i++){
        for(int j=0;j<COLS;j++){
            s[i][j]= 1;
        }
    }

    for(int ie=1;ie<80;ie++){
        temp[ie]=ie/10.0;
        arrayexp[0]=exp(-8/temp[ie]);
        arrayexp[1]=exp(8/temp[ie]);
        arrayexp[2]=exp(-4/temp[ie]);
        arrayexp[3]=exp(4/temp[ie]);
        arrayexp[4]=1.0;
        printf("\n%f",temp[ie]);
        // for(ii=0;ii<(N2*N2);ii++){
        //     x = rand()%ROWS;
        //     y = rand()%COLS;
        //     e1=energy(x,y,s);
        //     yrand = ((float)rand())/RAND_MAX;
        //     if(prob(-2*e1,arrayexp)>yrand){
        //         s[x][y]*= -1;
        //     }
        // }
        for(int _=0;_<N2;_++){ //10*N2*N2
            x = rand()%ROWS;
            y = rand()%COLS;
            e1=energy(x,y,s);
            yrand = ((float)rand())/RAND_MAX;
            if(prob(-2*e1,arrayexp)>yrand){
                s[x][y]*= -1;
            }
        }
        for(int ii=0;ii<NN;ii++){
        for(int _=0;_<2*NN;_++){
            x = rand()%ROWS;
            y = rand()%COLS;
            e1=energy(x,y,s);
            yrand = ((float)rand())/RAND_MAX;
            if(prob(-2*e1,arrayexp)>yrand){
                s[x][y]*= -1;
            }
        }
        MS=M(s);
        M_main+=MS;
        M2_main+=MS*MS;
        EnergyTotal=energy_t(s);
        energy_main+=EnergyTotal;
        energy2_main+=EnergyTotal*EnergyTotal;
        }
        printf("\t%f",M_main/NN);
        X=(1/temp[ie])*(M2_main/NN-pow((M_main/NN),2));
        printf("\t%.8f",X);
        printf("\t%f",energy_main/NN);
        Cv=(1/temp[ie])*(1/temp[ie])*(energy2_main/NN-pow((energy_main/NN),2));
        printf("\t%f",Cv);
        M_main=0.0;
        M2_main=0.0;
        energy_main=0.0;
        energy2_main=0.0;
    }
    return 0;
}