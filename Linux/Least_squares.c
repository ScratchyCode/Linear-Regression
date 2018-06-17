// Coded by ScratchyCode
// Find the best line that interpolates the entered data (considering negligible the uncertainties on the abscissas).
// Compile with -lm
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "stat.h"

#define INC 0.0001

int main(){
    int i, N;
    char c, file[100];
    
    // check that gnuplot is present on the system
    int gnupl = control();
    if(gnupl == 1){
        printf("\nYou need gnuplot to graph the results.");
        printf("\nInstall it with: sudo apt-get install gnuplot\n\n");
        exit(2);
    }
    
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
    printf("\nm = %.3lf\tsigma(m) = %.3lf\nq = %.3lf\tsigma(q) = %.3lf",M,sigmaM,Q,sigmaQ);
    printf("\n\nBest sigma(Y) = %.3lf",sigmaY);
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
    
    // creating fit
    printf("\nPlotting fit...\n");
    FILE *data = fopen("data.dat","w");
    if(data == NULL){
        perror("\nError");
        exit(1);
    }
    
    // writing experimental datas
    for(i=0; i<N; i++){
        fprintf(data,"%lf %lf %lf\n",x[i],y[i],errors[i]);
    }
    
    fclose(data);
    
    // creating fit points
    fit(M,Q,x,N);
    
    free(x);
    free(y);
       
    return 0;
}
