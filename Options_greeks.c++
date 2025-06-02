#include<iostream>
#include<math.h>
#include <cmath>    // for exp(), sqrt(), fabs(), M_PI
#define _USE_MATH_DEFINES
using namespace std;

// Standard normal PDF: φ(x) = exp(-x²/2) / √(2π)
double pdf(double x) {
    return (1.0 / sqrt(2.0 * M_PI)) * exp(-0.5 * x * x);
}

// Standard normal CDF via polynomial approximation
double cdf(double x) {
    double abs_x = fabs(x);
    double A = 1.0 / (1.0 + 0.2316419 * abs_x);
    double A_sum = A * (0.319381530 +
                       A * (-0.356563782 +
                       A * (1.781477937 +
                       A * (-1.821255978 +
                       1.330274429 * A))));
    double pdf_x = pdf(abs_x);
    double approx = 1.0 - pdf_x * A_sum;
    return (x >= 0.0) ? approx : (1.0 - approx);
}


double function(double s0, double r, double y, double T, double k ){
    double m0 = s0;
    double m1 = r;                          // y =sigma
    double m2 = y;
    double m3 = T;
    double m4 = sqrt(m3);   
    double m5 = (m4) *(m2);              
    double m6 = ((log(m0/k) + (m1*m3))/m5  + m5/2);
    double m7 = m6 - m5;
    double m8 = cdf(m6);
    double m9 = cdf(m7);
    double v = ((m0*m3) - k*exp(-m1*m3)*m9);
    return v;
}

void  adjoint(double s0, double r, double y, double T, double k, double v_dot){

    double m0 = s0;
    double m1 = r;                         
    double m2 = y;
    double m3 = T;
    double m4 = sqrt(m3);   
    double m5 = (m4) *(m2);              
    double m6 = ((log(m0/k) + (m1*m3))/m5  + m5/2);
    double m7 = m6 - m5;
    double m8 = cdf(m6);
    double m9 = cdf(m7);
    double v = ((m0*m3) - k*exp(-m1*m3)*m9);
    double m9_dot = ((-k*exp(-m1*m3))*v_dot);
    double m8_dot = m0*v_dot;
    double m7_dot = m9_dot * pdf(m7);
    double m6_dot = m8_dot * (pdf(m6) ) + m7_dot;
    double m5_dot = (-1) * ((m6_dot * (((log(m0/k) + m1*m3) / (m5*m5)) + 1/2 )) + m7_dot);
    double m4_dot = m5_dot * m2;
    double m3_dot = (m4_dot * (1/(2*sqrt(m3))) + ((m1/m5)*m6_dot) + (v_dot*k*m9*m1*exp(-m1*m3)));
    double m2_dot = m5_dot * m4;
    double m1_dot = (m6_dot * (m3/m5) + k*m9*m3*exp(-m1*m3));
    double m0_dot = (m6_dot * (1/(m5 * m0)) + m8*v_dot);
    
    cout << "Delta : " << m0_dot << endl;
    cout <<  "Rho : "  << m1_dot << endl;
    cout << "vega : " << m2_dot << endl;
    cout << "Theta : " << m3_dot << endl ;

}

int main(){
    double s0;
    double r; 
    double y; 
    double T; 
    double k; 
    double v_dot;
    cout << "Enter the value of S0 : ";
    cin >> s0 ;
    cout << "Enter the value of K : ";
    cin >> k;
    cout << "Enter the value of T : ";
    cin >> T;
    cout << "Enter the value of y : ";
    cin >> y;
    cout << "Enter the value of r : ";
    cin >> r;
    cout << "Enter the value of v_dot : ";
    cin >> v_dot ;
    
    adjoint(s0, r, y, T, k, v_dot);
    return 0;
}