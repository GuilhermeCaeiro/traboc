#!/bin/bash

### Experiments execution scheduler
#
#databases used in experiments
#4039187 a280.xml
#115959  att48.xml
#169678  brazil58.xml
#9166    burma14.xml
#12572180 d493.xml
#27159599 u724.xml

#Check groups to execute
if [[ "$1" == "1" || "$1" == "ALL" ]]; 
then
    #group1: experiments to compare performance among strategies christofides, factor2aprox, both, neareastneighbor, farthest 
    julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 none lagrangean
    julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 none lagrangean
    julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 none lagrangean
    julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 none lagrangean
    julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 none complementary
    julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 none complementary
    julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 none complementary
    julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 none complementary
    julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 2opt lagrangean
    julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 2opt lagrangean
    julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 2opt lagrangean
    julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 2opt lagrangean
    julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 2opt complementary
    julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 2opt complementary
    julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 2opt complementary
    julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 2opt complementary
    julia jgvtspsolver.jl factor2approximation ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 none lagrangean
    julia jgvtspsolver.jl factor2approximation ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 none lagrangean
    julia jgvtspsolver.jl factor2approximation ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 none lagrangean
    julia jgvtspsolver.jl factor2approximation ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 none lagrangean
    julia jgvtspsolver.jl factor2approximation ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 none complementary
    julia jgvtspsolver.jl factor2approximation ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 none complementary
    julia jgvtspsolver.jl factor2approximation ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 none complementary
    julia jgvtspsolver.jl factor2approximation ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 none complementary
    julia jgvtspsolver.jl factor2approximation ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 2opt lagrangean
    julia jgvtspsolver.jl factor2approximation ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 2opt lagrangean
    julia jgvtspsolver.jl factor2approximation ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 2opt lagrangean
    julia jgvtspsolver.jl factor2approximation ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 2opt lagrangean
    julia jgvtspsolver.jl factor2approximation ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 2opt complementary
    julia jgvtspsolver.jl factor2approximation ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 2opt complementary
    julia jgvtspsolver.jl factor2approximation ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 2opt complementary
    julia jgvtspsolver.jl factor2approximation ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 2opt complementary
    julia jgvtspsolver.jl christofidesandfactor2 ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 none lagrangean
    julia jgvtspsolver.jl christofidesandfactor2 ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 none lagrangean
    julia jgvtspsolver.jl christofidesandfactor2 ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 none lagrangean
    julia jgvtspsolver.jl christofidesandfactor2 ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 none lagrangean
    julia jgvtspsolver.jl christofidesandfactor2 ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 none complementary
    julia jgvtspsolver.jl christofidesandfactor2 ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 none complementary
    julia jgvtspsolver.jl christofidesandfactor2 ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 none complementary
    julia jgvtspsolver.jl christofidesandfactor2 ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 none complementary
    julia jgvtspsolver.jl christofidesandfactor2 ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 2opt lagrangean
    julia jgvtspsolver.jl christofidesandfactor2 ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 2opt lagrangean
    julia jgvtspsolver.jl christofidesandfactor2 ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 2opt lagrangean
    julia jgvtspsolver.jl christofidesandfactor2 ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 2opt lagrangean
    julia jgvtspsolver.jl christofidesandfactor2 ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 2opt complementary
    julia jgvtspsolver.jl christofidesandfactor2 ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 2opt complementary
    julia jgvtspsolver.jl christofidesandfactor2 ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 2opt complementary
    julia jgvtspsolver.jl christofidesandfactor2 ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 2opt complementary
    julia jgvtspsolver.jl nearestneighbor ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 none lagrangean
    julia jgvtspsolver.jl nearestneighbor ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 none lagrangean
    julia jgvtspsolver.jl nearestneighbor ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 none lagrangean
    julia jgvtspsolver.jl nearestneighbor ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 none lagrangean
    julia jgvtspsolver.jl nearestneighbor ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 none complementary
    julia jgvtspsolver.jl nearestneighbor ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 none complementary
    julia jgvtspsolver.jl nearestneighbor ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 none complementary
    julia jgvtspsolver.jl nearestneighbor ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 none complementary
    julia jgvtspsolver.jl nearestneighbor ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 2opt lagrangean
    julia jgvtspsolver.jl nearestneighbor ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 2opt lagrangean
    julia jgvtspsolver.jl nearestneighbor ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 2opt lagrangean
    julia jgvtspsolver.jl nearestneighbor ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 2opt lagrangean
    julia jgvtspsolver.jl nearestneighbor ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 2opt complementary
    julia jgvtspsolver.jl nearestneighbor ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 2opt complementary
    julia jgvtspsolver.jl nearestneighbor ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 2opt complementary
    julia jgvtspsolver.jl nearestneighbor ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 2opt complementary
    julia jgvtspsolver.jl farthestinsertion ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 none lagrangean
    julia jgvtspsolver.jl farthestinsertion ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 none lagrangean
    julia jgvtspsolver.jl farthestinsertion ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 none lagrangean
    julia jgvtspsolver.jl farthestinsertion ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 none lagrangean
    julia jgvtspsolver.jl farthestinsertion ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 none complementary
    julia jgvtspsolver.jl farthestinsertion ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 none complementary
    julia jgvtspsolver.jl farthestinsertion ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 none complementary
    julia jgvtspsolver.jl farthestinsertion ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 none complementary
    julia jgvtspsolver.jl farthestinsertion ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 2opt lagrangean
    julia jgvtspsolver.jl farthestinsertion ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 2opt lagrangean
    julia jgvtspsolver.jl farthestinsertion ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 2opt lagrangean
    julia jgvtspsolver.jl farthestinsertion ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 2opt lagrangean
    julia jgvtspsolver.jl farthestinsertion ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 2opt complementary
    julia jgvtspsolver.jl farthestinsertion ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 2opt complementary
    julia jgvtspsolver.jl farthestinsertion ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 2opt complementary
    julia jgvtspsolver.jl farthestinsertion ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 2opt complementary
fi 

if [[ "$1" == "2" || "$1" == "ALL" ]]; 
then
    #group2: best configuration from group1: change mi_function parameter
    julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 2opt complementary
    julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 current lbdecay1 50 0.5 2opt complementary
    julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 best lbdecay1 50 0.5 2opt complementary
    julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 5pct lbdecay1 50 0.5 2opt complementary
    julia jgvtspsolver.jl nearestneighbor ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 2opt complementary
    julia jgvtspsolver.jl nearestneighbor ./input/burma14.xml 10000 0.00001 2 current lbdecay1 50 0.5 2opt complementary
    julia jgvtspsolver.jl nearestneighbor ./input/burma14.xml 10000 0.00001 2 best lbdecay1 50 0.5 2opt complementary
    julia jgvtspsolver.jl nearestneighbor ./input/burma14.xml 10000 0.00001 2 5pct lbdecay1 50 0.5 2opt complementary
    julia jgvtspsolver.jl farthestinsertion ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 2opt complementary
    julia jgvtspsolver.jl farthestinsertion ./input/burma14.xml 10000 0.00001 2 current lbdecay1 50 0.5 2opt complementary
    julia jgvtspsolver.jl farthestinsertion ./input/burma14.xml 10000 0.00001 2 best lbdecay1 50 0.5 2opt complementary
    julia jgvtspsolver.jl farthestinsertion ./input/burma14.xml 10000 0.00001 2 5pct lbdecay1 50 0.5 2opt complementary
fi

if [[ "$1" == "3" || "$1" == "ALL" ]]; 
then
    #group3: best configuration from group1 applied to larger datasets
    julia jgvtspsolver.jl christofides ./input/att48.xml 10000 0.00001 2 best lbdecay1 50 0.5 2opt complementary
    julia jgvtspsolver.jl christofides ./input/brazil58.xml 10000 0.00001 2 best lbdecay1 50 0.5 2opt complementary
    julia jgvtspsolver.jl christofides ./input/a280.xml 10000 0.00001 2 best lbdecay1 50 0.5 2opt complementary
    julia jgvtspsolver.jl christofides ./input/d493.xml 10000 0.00001 2 best lbdecay1 50 0.5 2opt complementary
    julia jgvtspsolver.jl christofides ./input/u724.xml 10000 0.00001 2 best lbdecay1 50 0.5 2opt complementary
fi

if [[ "$1" == "4" || "$1" == "ALL" ]]; 
then
    #group4: parameters exploration using Gurobi and GLPK
    julia jgvtspsolver.jl gurobi ./input/burma14.xml 10000 secs none
    julia jgvtspsolver.jl gurobi ./input/burma14.xml 10000 combs none
    julia jgvtspsolver.jl gurobi ./input/burma14.xml 10000 both none
fi

if [[ "$1" == "5" || "$1" == "ALL" ]]; 
then
    #group5: optmization comparison using Gurobi and Christofides
    julia jgvtspsolver.jl gurobi ./input/burma14.xml 10000 both none
    julia jgvtspsolver.jl gurobi ./input/burma14.xml 10000 both christofides,0.00001,2,5pct,lbdecay1,50,0.5,2opt,complementary
    julia jgvtspsolver.jl gurobi ./input/att48.xml 10000 both none
    julia jgvtspsolver.jl gurobi ./input/att48.xml 10000 both christofides,0.00001,2,5pct,lbdecay1,50,0.5,2opt,complementary
    julia jgvtspsolver.jl gurobi ./input/brazil58.xml 10000 both none
    julia jgvtspsolver.jl gurobi ./input/brazil58.xml 10000 both christofides,0.00001,2,5pct,lbdecay1,50,0.5,2opt,complementary
fi

if [[ "$1" == "6" || "$1" == "ALL" ]]; 
then
    #group6: optmization comparison using Gurobi and Christofides, bigger datasets
    julia jgvtspsolver.jl gurobi ./input/a280.xml 10000 both none
    julia jgvtspsolver.jl gurobi ./input/a280.xml 10000 both christofides,0.00001,2,5pct,lbdecay1,50,0.5,2opt,complementary
    julia jgvtspsolver.jl gurobi ./input/d493.xml 10000 both none
    julia jgvtspsolver.jl gurobi ./input/d493.xml 10000 both christofides,0.00001,2,5pct,lbdecay1,50,0.5,2opt,complementary
fi    