// Coded by ScratchyCode
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define INC 0.0001

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

void extrapolation(int N, double *x, double *errors, double sigmaY, double M, double Q){
    int i;
    double alpha, pointX, pointY, sigmaPointY, xMean=0, sumXnum=0, sumXdenom=0;
    
    printf("Insert the point's abscissa: ");
    scanf("%lf",&pointX);
    
    // xMean
    for(i=0; i<N; i++){
        // xMean
        sumXnum += x[i] / pow(errors[i],2);
        sumXdenom += 1 / pow(errors[i],2);
    }
    
    xMean = sumXnum/sumXdenom;
    
    pointY = (M * pointX) + Q;
    alpha = sqrt(sigma(x,N,mean(x,N)) + (N/pow(sigmaY,2)));
    sigmaPointY = (1/alpha) * sqrt(pow(pointX,2) - (2*pointX*xMean) + pow(pointX,2));
    
    printf("F(%.10lf) = %.10lf +- %.10lf\n",pointX,pointY,sigmaPointY);
    
    return ;
}

double Min(double array[], int dim){
    long long int i;
    double min=1E20;
    
    for(i=0; i<dim; i++){
        if(min > array[i]){
            min = array[i];
        }
    }
    
    return min;
}

double Max(double array[], int dim){
    long long int i;
    double max=-1E20;
    
    for(i=0; i<dim; i++){
        if(max < array[i]){
            max = array[i];
        }
    }
    
    return max;
}

void fit(double M, double Q, double array[], int dim){
    double x, y, inf, sup;
    FILE *fit = fopen("fit.dat","w");
    if(fit == NULL){
        perror("\nError");
        exit(1);
    }
    
    inf = Min(array,dim);
    sup = Max(array,dim);
    x = inf - 0.5;
    
    do{
        y = M * x + Q;
        fprintf(fit,"%lf %lf\n",x,y);
        x += INC;
    }while(x <= sup + 0.5);
    
    fflush(fit);
    fclose(fit);
    
    return ;
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

int control(void){
    char path[] = "/usr/bin/gnuplot";
    
    FILE *pf = fopen(path,"r");
    if(pf == NULL){
        fclose(pf);
        return 1;
    }else{
        fclose(pf);
        return 0;
    }
}
