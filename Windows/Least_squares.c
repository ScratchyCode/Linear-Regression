// Coded by ScratchyCode
// Find the best line that interpolates the entered data (considering negligible the uncertainties on the abscissas).
// Compile with -lm
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define INC 0.0001

double Mbest(int N, double *x, double *y);
double Qbest(int N, double *x, double *y);
double bestSigma(int N, double *x, double *y, double M, double Q);
double uM(int N, double sigmaY, double *x);
double uQ(int N, double sigmaY, double *x);
double covariance(double *x, double *y, int N);
double correlation(double *x, double *y, int N);
double mean(double array[], int N);
double sigma(double array[], int N, double mean);
int linesFile(char file[]);

int main(){
    int i, N;
    double M, Q, error;
    char c, file[100];
    
    printf("\nEnter the name or path of data's file: ");
    fgets(file,sizeof(file),stdin);
    file[strlen(file)-1] = '\0';
    
    printf("Enter the uncertainty regarding the entered data: ");
    scanf("%lf",&error);
    
    // the number of lines in the file is the number of datas to be saved
    N = linesFile(file);
    if(N <= 2){
        printf("\nError: insufficient data number.\n");
        exit(2);
    }
    
    // creating the data's arrays
    double *x = calloc(N,sizeof(double));
    double *y = calloc(N,sizeof(double));
    if(x == NULL || y == NULL){
            perror("\nerror");
            printf("\n");
            exit(1);
    }
    
    FILE *inputFile = fopen(file,"r");
    if(inputFile == NULL){
        perror("\nError");
        exit(1);
    }
    
    // reading from file
    for(i=0; i<N; i++){
        fscanf(inputFile,"%lf %lf\n",&x[i],&y[i]);
    }
    
    fclose(inputFile);
    
    // determine linear coefficients
    M = Mbest(N,x,y);
    Q = Qbest(N,x,y);
    
    // defining best sigma(Y) and correlation coefficient
    double sigmaY = bestSigma(N,x,y,M,Q);
    double sigmaM = uM(N,sigmaY,x);
    double sigmaQ = uQ(N,sigmaY,x);
    double cov = covariance(x,y,N);
    double cor = correlation(x,y,N);
    
    //printf("\nCov(X,Y) = %lf",cov);
    //printf("\nCor(X,Y) = %.3lf",cor);
    printf("\nBest sigma(Y) = %.3lf\n",sigmaY);
    printf("\nThe best line Y = mX + q that fit data is:");
    printf("\nm = %.3lf\tsigma(m) = %.3lf\nq = %.3lf\tsigma(q) = %.3lf\n\n",M,Q,sigmaM,sigmaQ);
       
    return 0;
}

double Mbest(int N, double *x, double *y){
    int i;
    double Msum1=0, Msum2=0, Msum3=0, Msum4=0;
    
    for(i=0; i<N; i++){
        Msum1 = Msum1 + x[i] * y[i];
        Msum2 = Msum2 + x[i];
        Msum3 = Msum3 + y[i];
        Msum4 = Msum4 + pow(x[i],2);
    }
    
    return ((N*Msum1 - Msum2*Msum3) / (N*Msum4 - pow(Msum2,2)));   
}

double Qbest(int N, double *x, double *y){
    int i;
    double Qsum1=0, Qsum2=0, Qsum3=0, Qsum4=0;
    
    for(i=0; i<N; i++){
        Qsum1 = Qsum1 + y[i];
        Qsum2 = Qsum2 + pow(x[i],2);
        Qsum3 = Qsum3 + x[i];
        Qsum4 = Qsum4 + x[i] * y[i];
    }
    
    return ((Qsum1*Qsum2 - Qsum3*Qsum4) / (N*Qsum2 - pow(Qsum3,2)));
}

double bestSigma(int N, double *x, double *y, double M, double Q){
    int i;
    double sum=0;
    
    for(i=0; i<N; i++){
        sum += pow(y[i] - Q - M * x[i],2);
    }
    
    return sqrt(sum/(N-2));
}

double uM(int N, double sigmaY, double *x){
    int i;
    double sum1=0, sum2=0, delta;
    
    for(i=0; i<N; i++){
        sum1 += pow(x[i],2);
        sum2 += x[i];
    }
    
    sum2 = pow(sum2,2);
    delta = N * sum1 - sum2;
    
    return sigmaY * sqrt((sum1/delta));
}

double uQ(int N, double sigmaY, double *x){
    int i;
    double sum1=0, sum2=0, delta;
    
    for(i=0; i<N; i++){
        sum1 += pow(x[i],2);
        sum2 += x[i];
    }
    
    sum2 = pow(sum2,2);
    delta = (N * sum1) - sum2;
    
    return sigmaY * sqrt(N/delta);
}

double covariance(double *x, double *y, int N){
    int i;
    double meanX = mean(x,N);
    double meanY = mean(y,N);
    double sum=0;
    
    for(i=0; i<N; i++){
        sum += (x[i] - meanX) * (y[i] - meanY);
    }
    
    return sum/N;
}

double correlation(double *x, double *y, int N){
    double meanX = mean(x,N);
    double meanY = mean(y,N);
    double sigmaX = sigma(x,N,meanX);
    double sigmaY = sigma(y,N,meanY);
    double cov = covariance(x,y,N);
    
    return cov/(sigmaX * sigmaY);
}

double mean(double array[], int N){
    int i;
    double sum=0;
    
    for(i=0; i<N; i++){
        sum += array[i];
    }
    
    return sum/N;
}

double sigma(double array[], int N, double mean){
    int i;
    double sigma, sum=0;
    
    for(i=0; i<N; i++){
        sum += pow((array[i] - mean),2);
    }
    
    sigma = sqrt(sum/(N-1));
    
    return sigma;
}

int linesFile(char file[]){
    int lines=0;
    char c;
    
    FILE *input = fopen(file,"r");
    if(input == NULL){
        perror("\nError");
        exit(1);
    }
    
    while((c = getc(input)) != EOF){
        if(c == '\n'){
            lines++;
        }
    }
    
    fclose(input);
    
    return lines;
}
