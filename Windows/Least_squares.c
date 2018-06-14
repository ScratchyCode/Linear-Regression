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
double linearParamCovariance(double meanX, double sigmaX, double sigmaY, int N, double M, double Q);
double linearParamCorrelation(double *x, int N);
double mean(double array[], int N);
double sigma(double array[], int N, double mean);
double sign(double x);
int linesFile(char file[]);

int main(){
    int i, N;
    double M, Q, error;
    char c, file[100];
    
    printf("Enter the name or path of file: ");
    fgets(file,sizeof(file),stdin);
    file[strlen(file)-1] = '\0';
    
    //printf("Enter the experimental uncertainty sigma(Y): ");
    //scanf("%lf",&error);
    
    // the file's lines number is the number of points to be saved
    N = linesFile(file);
    if(N <= 2){
        printf("\nError: insufficient data number.\n");
        exit(2);
    }
    
    // creating data's arrays
    double *x = calloc(N,sizeof(double));
    double *y = calloc(N,sizeof(double));
    if(x == NULL || y == NULL){
            perror("\nerror");
            printf("\n");
            exit(1);
    }
    
    // reading from file
    FILE *inputFile = fopen(file,"r");
    if(inputFile == NULL){
        perror("\nError");
        exit(1);
    }
    
    for(i=0; i<N; i++){
        fscanf(inputFile,"%lf %lf\n",&x[i],&y[i]);
    }
    
    fclose(inputFile);
    
    // determine linear coefficients
    M = Mbest(N,x,y);
    Q = Qbest(N,x,y);
    
    // defining best sigma(Y) and correlation coefficient
    double sigmaY = fabs(bestSigma(N,x,y,M,Q)); // <-- residuals analysis
    double sigmaM = uM(N,sigmaY,x);
    double sigmaQ = uQ(N,sigmaY,x);
    double cov = covariance(x,y,N);
    double cor = correlation(x,y,N);
    double lCov = linearParamCovariance(mean(x,N),sigma(x,N,mean(x,N)),sigmaY,N,M,Q);
    double lCor = linearParamCorrelation(x,N);
    
    printf("\nThe best linear fit Y = mX + q is:");
    printf("\nm = %.3lf\tsigma(m) = %.3lf\nq = %.3lf\tsigma(q) = %.3lf",M,Q,sigmaM,sigmaQ);
    printf("\nBest sigma(Y) = %.3lf\n",sigmaY);
    printf("\nCov(X,Y) = %.3lf",cov);
    printf("\nCor(X,Y) = %.3lf",cor);
    printf("\nCov(m,c) = %.3lf",lCov);
    printf("\nCor(m,c) = %.3lf",lCor);
    
    // interpolation and extrapolation
    int choice;
    double pointX, pointY, sigmaPointY, alpha;
    printf("\n\nDo you want to extrapolate a point with the calculated linear regression? (1 = YES | 0 = NO): ");
    scanf("%d",&choice);
    if(choice != 1){
        exit(0);
    }
    
    printf("Insert the point's abscissa: ");
    scanf("%lf",&pointX);
    
    pointY = (M * pointX) + Q;
    alpha = sqrt(sigma(x,N,mean(x,N)) + (N/pow(sigmaY,2)));
    sigmaPointY = (1/alpha) * sqrt(pow(pointX,2) - (2*pointX*mean(x,N)) + pow(pointX,2));
    
    printf("F(%.3lf) = %.3lf +- %.3lf\n",pointX,pointY,sigmaPointY);
       
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
        sum += ((x[i] - meanX) * (y[i] - meanY));
    }
    
    return sum/(N-1);
}

double correlation(double *x, double *y, int N){
    double meanX = mean(x,N);
    double meanY = mean(y,N);
    double sigmaX = sigma(x,N,meanX);
    double sigmaY = sigma(y,N,meanY);
    double cov = covariance(x,y,N);
    
    return cov/(sigmaX * sigmaY);
}

double linearParamCovariance(double meanX, double sigmaX, double sigmaY, int N, double M, double Q){
    return -( ((meanX)/(pow(sigmaX,2))) * (pow(sigmaY,2)/N) );
}

double linearParamCorrelation(double *x, int N){
    double meanX = mean(x,N);
    double sigmaX = sigma(x,N,meanX);
    
    return -( sign(meanX) / sqrt(1 + (pow(sigmaX,2)/meanX)) ) ;
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

double sign(double x){
    if(x >= 0){
        return 1;
    }else if(x < 0){
        return -1;
    }
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
