// Coded by ScratchyCode

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

void fit(double M, double Q, double array[], int dim);

double Min(double array[], int dim);

double Max(double array[], int dim);

double sign(double x);

int linesFile(char file[]);

int control(void);
