#!/bin/bash
#databases
#4039187 a280.xml
#115959  att48.xml
#169678  brazil58.xml
#9166    burma14.xml
#12572180 d493.xml
#27159599 u724.xml

#
#
#step1: experiments to compare performance among strategies christofides, factor2aprox, both, neareastneighbor, farthest 
julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 none lagrangean
julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 none lagrangean
julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 none lagrangean
julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 none lagrangean
julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 none complementary
julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 none complementary
julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 none complementary
julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 none complementary
julia jgvtspsolver.jl factor2approximation ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 none lagrangean
julia jgvtspsolver.jl factor2approximation ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 none lagrangean
julia jgvtspsolver.jl factor2approximation ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 none lagrangean
julia jgvtspsolver.jl factor2approximation ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 none lagrangean
julia jgvtspsolver.jl factor2approximation ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 none complementary
julia jgvtspsolver.jl factor2approximation ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 none complementary
julia jgvtspsolver.jl factor2approximation ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 none complementary
julia jgvtspsolver.jl factor2approximation ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 none complementary
julia jgvtspsolver.jl christofidesandfactor2 ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 none lagrangean
julia jgvtspsolver.jl christofidesandfactor2 ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 none lagrangean
julia jgvtspsolver.jl christofidesandfactor2 ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 none lagrangean
julia jgvtspsolver.jl christofidesandfactor2 ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 none lagrangean
julia jgvtspsolver.jl christofidesandfactor2 ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 none complementary
julia jgvtspsolver.jl christofidesandfactor2 ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 none complementary
julia jgvtspsolver.jl christofidesandfactor2 ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 none complementary
julia jgvtspsolver.jl christofidesandfactor2 ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 none complementary
julia jgvtspsolver.jl nearestneighbor ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 none lagrangean
julia jgvtspsolver.jl nearestneighbor ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 none lagrangean
julia jgvtspsolver.jl nearestneighbor ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 none lagrangean
julia jgvtspsolver.jl nearestneighbor ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 none lagrangean
julia jgvtspsolver.jl nearestneighbor ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 none complementary
julia jgvtspsolver.jl nearestneighbor ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 none complementary
julia jgvtspsolver.jl nearestneighbor ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 none complementary
julia jgvtspsolver.jl nearestneighbor ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 none complementary
julia jgvtspsolver.jl farthestinsertion ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 none lagrangean
julia jgvtspsolver.jl farthestinsertion ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 none lagrangean
julia jgvtspsolver.jl farthestinsertion ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 none lagrangean
julia jgvtspsolver.jl farthestinsertion ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 none lagrangean
julia jgvtspsolver.jl farthestinsertion ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 none complementary
julia jgvtspsolver.jl farthestinsertion ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 none complementary
julia jgvtspsolver.jl farthestinsertion ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 none complementary
julia jgvtspsolver.jl farthestinsertion ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 none complementary

#step2: best configuration from step1: change mi_function parameter
#julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 current lbdecay2 50 0.5 none lagrangean
#julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 best lbdecay2 50 0.5 none lagrangean
#julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 5pct lbdecay2 50 0.5 none lagrangean

#step3: best configuration from step1 applied to larger datasets
#julia jgvtspsolver.jl christofides ./input/att48.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 none lagrangean
#julia jgvtspsolver.jl christofides ./input/brazil58.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 none lagrangean
#julia jgvtspsolver.jl christofides ./input/a280.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 none lagrangean
#julia jgvtspsolver.jl christofides ./input/d493.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 none lagrangean
