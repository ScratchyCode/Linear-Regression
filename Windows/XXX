// Coded by ScratchyCode
// Find the best line that interpolates the entered data (considering negligible the uncertainties on the abscissas).
// Compile with -lm
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define INC 0.0001

double mean(double array[], int N);
double sigma(double array[], int N, double mean);
double covariance(int N, double *x, double *y);
double correlation(int N, double *x, double *y);
double linearParamCovariance(int N, double meanX, double sigmaX, double sigmaY, double M, double Q);
double linearParamCorrelation(int N, double *x);
double Mbest(int N, double *x, double *y, double *errors);
double Qbest(int N, double *x, double *y, double *errors, double M);
double uM(int N, double *x, double *y, double *errors, double M, double Q);
double uQ(int N, double *x, double *errors, double sigmaM);
double bestSigma(int N, double *x, double *y, double M, double Q);
void extrapolation(int N, double *x, double sigmaY, double M, double Q);
int linesFile(char file[]);

int main(){
    int i, N;
    char c, file[100];
    
    printf("Enter the name or path of file: ");
    fgets(file,sizeof(file),stdin);
    file[strlen(file)-1] = '\0';
    
    // the file's lines number is the number of points to be saved
    N = linesFile(file);
    if(N <= 2){
        printf("\nError: insufficient data number.\n");
        exit(2);
    }
    
    // creating data's arrays
    double *x = calloc(N,sizeof(double));
    double *y = calloc(N,sizeof(double));
    double *errors = calloc(N,sizeof(double));
    if(x == NULL || y == NULL || errors == NULL){
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
        fscanf(inputFile,"%lf %lf %lf\n",&x[i],&y[i],&errors[i]);
    }
    
    fclose(inputFile);
    
    // determine linear coefficients
    double M = Mbest(N,x,y,errors);
    double Q = Qbest(N,x,y,errors,M);
    double sigmaM = fabs(uM(N,x,y,errors,M,Q));
    double sigmaQ = fabs(uQ(N,x,errors,sigmaM));
    
    // defining best sigma(Y) and correlation coefficient
    double sigmaY = fabs(bestSigma(N,x,y,M,Q)); // <-- residuals analysis
    double cov = covariance(N,x,y);
    double cor = correlation(N,x,y);
    double lCov = linearParamCovariance(N,mean(x,N),sigma(x,N,mean(x,N)),sigmaY,M,Q);
    double lCor = linearParamCorrelation(N,x);
    
    // Chi-square test
    int freedomDegrees = N - 2; // infer 2 parameters (): M and Q
    double chi2 = 0, rChi2;
    
    for(i=0; i<N; i++){
        chi2 += pow(y[i] - ((M * x[i]) + Q),2) / pow(errors[i],2);
    }
    
    rChi2 = chi2 / freedomDegrees;
    
    printf("\nThe best linear fit Y = mX + q is:");
    printf("\nm = %.3lf\tsigma(m) = %.3lf\nq = %.3lf\tsigma(q) = %.3lf",M,Q,sigmaM,sigmaQ);
    printf("\nBest sigma(Y) = %.3lf\n",sigmaY);
    printf("\nCov(X,Y) = %.3lf",cov);
    printf("\nCor(X,Y) = %.3lf",cor);
    printf("\nCov(m,c) = %.3lf",lCov);
    printf("\nCor(m,c) = %.3lf",lCor);
    printf("\nChi square = %.3lf",chi2);
    printf("\nReduced Chi square = %.3lf",rChi2);
    
    // interpolation and extrapolation
    int choice;
    double pointX, pointY, sigmaPointY, alpha;
    printf("\n\nDo you want to extrapolate a point with the calculated linear regression? (1 = YES | 0 = NO): ");
    scanf("%d",&choice);
    if(choice == 1){
        extrapolation(N,x,sigmaY,M,Q);
    }
    
    free(x);
    free(y);
       
    return 0;
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

double covariance(int N, double *x, double *y){
    int i;
    double meanX = mean(x,N);
    double meanY = mean(y,N);
    double sum=0;
    
    for(i=0; i<N; i++){
        sum += ((x[i] - meanX) * (y[i] - meanY));
    }
    
    return sum/(N-1);
}

double correlation(int N, double *x, double *y){
    double meanX = mean(x,N);
    double meanY = mean(y,N);
    double sigmaX = sigma(x,N,meanX);
    double sigmaY = sigma(y,N,meanY);
    double cov = covariance(N,x,y);
    
    return cov/(sigmaX * sigmaY);
}

double linearParamCovariance(int N, double meanX, double sigmaX, double sigmaY, double M, double Q){
    return -( ((meanX)/(pow(sigmaX,2))) * (pow(sigmaY,2)/N) );
}

double linearParamCorrelation(int N, double *x){
    double meanX = mean(x,N);
    double sigmaX = sigma(x,N,meanX);
    
    return -(meanX / sqrt(pow(sigmaX,2) + pow(meanX,2)));
}

double Mbest(int N, double *x, double *y, double *errors){
    int i;
    double xMean=0, yMean=0, x2Mean=0, xyMean=0;
    double sumXnum=0, sumYnum=0, sumX2num=0, sumXYnum=0;
    double sumXdenom=0, sumYdenom=0, sumX2denom=0, sumXYdenom=0;
    
    for(i=0; i<N; i++){
        // xMean
        sumXnum += x[i] / pow(errors[i],2);
        sumXdenom += 1 / pow(errors[i],2);
        // yMean
        sumYnum += y[i] / pow(errors[i],2);
        sumYdenom += 1 / pow(errors[i],2);
        // x2Mean
        sumX2num += pow(x[i],2) / pow(errors[i],2);
        sumX2denom += 1 / pow(errors[i],2);
        // xyMean
        sumXYnum += (x[i] * y[i]) / pow(errors[i],2);
        sumXYdenom += 1 / pow(errors[i],2);
    }
    
    xMean = sumXnum / sumXdenom;
    yMean = sumYnum / sumYdenom;
    x2Mean = sumX2num / sumX2denom;
    xyMean = sumXYnum / sumXYdenom;
    
    return ((xyMean - (xMean * yMean)) / (x2Mean - pow(xMean,2)));
}

double Qbest(int N, double *x, double *y, double *errors, double M){
    int i;
    double xMean=0, yMean=0;
    double sumXnum=0, sumYnum=0;
    double sumXdenom=0, sumYdenom=0;
    
    for(i=0; i<N; i++){
        // xMean
        sumXnum += x[i] / pow(errors[i],2);
        sumXdenom += 1 / pow(errors[i],2);
        // yMean
        sumYnum += y[i] / pow(errors[i],2);
        sumYdenom += 1 / pow(errors[i],2);
    }
    
    xMean = sumXnum / sumXdenom;
    yMean = sumYnum / sumYdenom;
    
    return (yMean - (M * xMean));
}

double uM(int N, double *x, double *y, double *errors, double M, double Q){
    int i;
    double var=0, sumXnum=0, sumXdenom=0, sumX2num=0,sumX2denom=0, xMean=0, x2Mean=0;
    
    // var
    for(i=0; i<N; i++){
        var += pow(y[i] - ((M * x[i]) + Q),2) / pow(errors[i],2);
    }
    
    var /= (N - 2);
    
    // xMean and x2Mean
    for(i=0; i<N; i++){
        // xMean
        sumXnum += x[i] / pow(errors[i],2);
        sumXdenom += 1 / pow(errors[i],2);
        // x2Mean
        sumX2num += pow(x[i],2) / pow(errors[i],2);
        sumX2denom += 1 / pow(errors[i],2);
    }
    
    xMean = sumXnum/sumXdenom;
    x2Mean = sumX2num/sumX2denom;
    
    return sqrt(var / (N * (x2Mean - pow(xMean,2))));
}

double uQ(int N, double *x, double *errors, double sigmaM){
    int i;
    double sumX2num=0, sumX2denom=0, x2Mean=0;
    
    //x2Mean
    for(i=0; i<N; i++){
        // x2Mean
        sumX2num += pow(x[i],2) / pow(errors[i],2);
        sumX2denom += 1 / pow(errors[i],2);
    }
    
    x2Mean = sumX2num / sumX2denom;
    
    return sqrt(pow(sigmaM,2) * x2Mean);
}

double bestSigma(int N, double *x, double *y, double M, double Q){
    int i;
    double sum=0;
    
    for(i=0; i<N; i++){
        sum += pow(y[i] - Q - M * x[i],2);
    }
    
    return sqrt(sum/(N-2));
}

void extrapolation(int N, double *x, double sigmaY, double M, double Q){
    double alpha, pointX, pointY, sigmaPointY;
    
    printf("Insert the point's abscissa: ");
    scanf("%lf",&pointX);
    
    pointY = (M * pointX) + Q;
    alpha = sqrt(sigma(x,N,mean(x,N)) + (N/pow(sigmaY,2)));
    sigmaPointY = (1/alpha) * sqrt(pow(pointX,2) - (2*pointX*mean(x,N)) + pow(pointX,2));
    
    printf("F(%.3lf) = %.3lf +- %.3lf\n",pointX,pointY,sigmaPointY);
    
    return ;
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
