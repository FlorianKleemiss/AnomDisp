#pragma once
#include <cmath>
#include <vector>
#include "convenience.h"

double alpha_coef(int l, int m, int mp, double t0, double alpha){
    if (m>l || mp > l) return 0.0;
    double st0 = sin(t0/2.0);
    double ct0 = cos(t0/2.0);
    if (mp == 0){
        double part1 = cos(m*alpha) * ft[l] / ft[l-m] * pow(ct0,2*l);
        double sum = 0;
        for(int rho=0; rho <= l-m; rho++)
            sum += pow(-1,rho) * binom(l-m,rho) * binom(l+m,l-rho) * pow(st0/ct0,m+2*rho);
        return part1 * sum;
    }
    double part1 = cos(m*alpha) * ft[l-mp] / ft[l-m];
    if (part1 == 0) return 0.0;
    double sum = 0;
    for(int rho=0; rho <= l-m; rho++){
      int f0 = binom(l-m,rho);
      int f1 = binom(l+m,m+mp+rho);
      int f2 = binom(l+m,m-mp+rho);
      double q1, q2, s;
      if (f1 != 0 && f2 != 0){
        q1 = pow(st0,m+mp+2*rho) * pow(ct0,2*l-m-2*rho-mp);
        q2 = pow(st0,m-mp+2*rho) * pow(ct0,2*l-m-2*rho+mp);
        if (q1 != 0 && q2 != 0)
          s = (pow(-1,mp) * f1 * q1 + f2 * q2);
        else if (q1 != 0)
          s = (pow(-1,mp) * f1 * q1);
        else if (q2 != 0)
          s = f2*q2;
        else
          s = 0;
      }
      else if (f1 == 0){
        double temp = m-mp+2*rho;
        if (temp < 0 && st0 == 0)
          q2 = 0;
        else
          q2 = pow(st0,m-mp+2*rho) * pow(ct0,2*l-m-2*rho+mp);
        s = f2 * q2;
      }
      else if (f2 == 0){
        q1 = pow(st0,m+mp+2*rho) * pow(ct0,2*l-m-2*rho-mp);
        s = pow(-1,mp) * f1 * q1;
      }
      else
        s = 0;
      sum += pow(-1,rho) * f0 * s;
    }
    return part1 * sum;
}

double alpha_bar_coef(int l, int m, int mp, double t0, double alpha){
    if (l-m < 0) return 0.0;
    double st0 = sin(t0/2.0);
    double ct0 = cos(t0/2.0);
    if (mp == 0){
        double part1 = sin(m*alpha) * ft[l] / ft[l-m] * pow(ct0,2*l);
        double sum = 0;
        for(int rho=0; rho <= l-m; rho++)
            sum += pow(-1,rho) * binom(l-m,rho) * binom(l+m,l-rho) * pow(st0/ct0,m+2*rho);
        return part1 * sum;
    }
    double part1 = sin(m*alpha) * ft[l-mp] / ft[l-m];
    if (part1 == 0) return 0.0;
    double sum = 0;
    for(int rho=0; rho <= l-m; rho++){
      int f0 = binom(l-m,rho);
      int f1 = binom(l+m,m+mp+rho);
      int f2 = binom(l+m,m-mp+rho);
      double q1, q2, s;
      if (f1 != 0 && f2 != 0){
        q1 = pow(st0,m+mp+2*rho) * pow(ct0,2*l-m-2*rho-mp);
        q2 = pow(st0,m-mp+2*rho) * pow(ct0,2*l-m-2*rho+mp);
        if (q1 != 0 && q2 != 0)
          s = (pow(-1,mp) * f1 * q1 + f2 * q2);
        else if (q1 != 0)
          s = (pow(-1,mp) * f1 * q1);
        else if (q2 != 0)
          s = f2*q2;
        else
          s = 0;
      }
      else if (f1 == 0){
        q2 = pow(st0,m-mp+2*rho) * pow(ct0,2*l-m-2*rho+mp);
        s = f2 * q2;
      }
      else if (f2 == 0){
        q1 = pow(st0,m+mp+2*rho) * pow(ct0,2*l-m-2*rho-mp);
        s = pow(-1,mp) * f1 * q1;
      }
      else
        s = 0;
      sum += pow(-1,rho) * f0 * s;
    }
    return part1 * sum;
}

double beta_coef(int l, int m, int mp, double t0, double alpha){
    if (mp == 0) return 0.0;
    if (m == 0) return 0.0;
    if (l-m < 0) return 0.0;
    double st0 = sin(t0/2.0);
    double ct0 = cos(t0/2.0);
    double part1 = sin(m*alpha) * ft[l-mp] / ft[l-m];
    if (part1 == 0) return 0.0;
    double sum = 0;
    for(int rho=0; rho <= l-m; rho++){
      int f0 = binom(l-m,rho);
      int f1 = binom(l+m,m+mp+rho);
      int f2 = binom(l+m,m-mp+rho);
      double q1, q2, s;
      if (f1 != 0 && f2 != 0){
        q1 = pow(st0,m+mp+2*rho) * pow(ct0,2*l-m-2*rho-mp);
        q2 = pow(st0,m-mp+2*rho) * pow(ct0,2*l-m-2*rho+mp);
        if (q1 != 0 && q2 != 0)
          s = (pow(-1,mp) * f1 * q1 - f2 * q2);
        else if (q1 != 0)
          s = (pow(-1,mp) * f1 * q1);
        else if (q2 != 0)
          s = -f2*q2;
        else
          s = 0;
      }
      else if (f1 == 0){
        q2 = pow(st0,m-mp+2*rho) * pow(ct0,2*l-m-2*rho+mp);
        s = -f2 * q2;
      }
      else if (f2 == 0){
        q1 = pow(st0,m+mp+2*rho) * pow(ct0,2*l-m-2*rho-mp);
        s = pow(-1,mp) * f1 * q1;
      }
      else
        s = 0;
      sum += pow(-1,rho) * f0 * s;
    }
    return part1 * sum;
}

double beta_bar_coef(int l, int m, int mp, double t0, double alpha){
    if (mp == 0) return 0.0;
    if (l-m < 0) return 0.0;
    double st0 = sin(t0/2.0);
    double ct0 = cos(t0/2.0);
    double part1 = cos(m*alpha) * ft[l-mp] / ft[l-m];
    if (part1 == 0) return 0.0;
    double sum = 0;
    for(int rho=0; rho <= l-m; rho++){
      int f0 = binom(l-m,rho);
      int f1 = binom(l+m,m+mp+rho);
      int f2 = binom(l+m,m-mp+rho);
      double q1, q2, s;
      if (f1 != 0 && f2 != 0){
        q1 = pow(st0,m+mp+2*rho) * pow(ct0,2*l-m-2*rho-mp);
        q2 = pow(st0,m-mp+2*rho) * pow(ct0,2*l-m-2*rho+mp);
        if (q1 != 0 && q2 != 0)
          s = (pow(-1,mp+1) * f1 * q1 + f2 * q2);
        else if (q1 != 0)
          s = (pow(-1,mp+1) * f1 * q1);
        else if (q2 != 0)
          s = -f2*q2;
        else
          s = 0;
      }
      else if (f1 == 0){
        q2 = pow(st0,m-mp+2*rho) * pow(ct0,2*l-m-2*rho+mp);
        s = -f2 * q2;
      }
      else if (f2 == 0){
        q1 = pow(st0,m+mp+2*rho) * pow(ct0,2*l-m-2*rho-mp);
        s = pow(-1,mp) * f1 * q1;
      }
      else
        s = 0;
      sum += pow(-1,rho) * f0 * s;
    }
    return part1 * sum;
}