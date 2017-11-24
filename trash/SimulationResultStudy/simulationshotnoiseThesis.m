clear all
close all
clc

load physconst

PdBm=-75:1:0;
PW=10.^((PdBm-30)/10);
B=30e9;
lambda=1550e-9;

SNU=sqrt(PW*B*h*c/lambda)*10^6