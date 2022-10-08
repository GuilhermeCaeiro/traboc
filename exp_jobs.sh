#!/bin/bash

### Experiments execution scheduler
#
#methods used in experiments:
#christofides factor2approximation nearestneighbor farthestinsertion
#
#databases used in experiments
#4039187 a280.xml
#115959  att48.xml
#169678  brazil58.xml
#822029  bier127.xml
#9166    burma14.xml
#12572180 d493.xml
#27159599 u724.xml
#
#mi_functions
#current best 5pct 1pct
#
#epsilon_strategies
#lbdecay1 lbdecay2 itdecay
#
#epsilon_decay_intervals
#30 50
#
#epsilon_decay_multipliers
#0.25 0.5 0.75
#
#local_search_algorithms
#2opt none
#
#primal_inputs
#lagrangean complementary

#Check groups to execute
if [[ "$1" == "L1" || "$1" == "ALL" ]]; 
then
    #group1: experiments to compare performance among strategies christofides, factor2aprox, both, neareastneighbor, farthest 
    methods="christofides factor2approximation nearestneighbor farthestinsertion"
    databases="burma14.xml att48.xml brazil58.xml bier127.xml"
    mi_functions="current best 5pct 1pct"
    epsilon_strategies="lbdecay1 lbdecay2 itdecay"
    epsilon_decay_intervals="30 50"
    epsilon_decay_multipliers="0.25 0.5 0.75"
    local_search_algorithms="2opt none"
    primal_inputs="lagrangean complementary"
    #
    for method in $methods; 
    do
        for database in $databases;
        do
            for mi_function in $mi_functions;
            do
                for epsilon_strategy in $epsilon_strategies;
                do
                    if [[ "$epsilon_strategy" == "static" ]];
                    then
                        for local_search_algorithm in $local_search_algorithms;
                        do
                            for primal_input in $primal_inputs;
                            do
                                #exec experiment
                                julia jgvtspsolver.jl $method ./input/$database 10000 0.00001 2 $mi_function $epsilon_strategy 0 0 $local_search_algorithm $primal_input
                            done
                        done
                    else
                        for epsilon_decay_interval in $epsilon_decay_intervals;
                        do
                            for epsilon_decay_multiplier in $epsilon_decay_multipliers;
                            do
                                for local_search_algorithm in $local_search_algorithms; 
                                do
                                    for primal_input in $primal_inputs; 
                                    do
                                        #exec experiment
                                        julia jgvtspsolver.jl $method ./input/$database 10000 0.00001 2 $mi_function $epsilon_strategy $epsilon_decay_interval $epsilon_decay_multiplier $local_search_algorithm $primal_input  
                                    done
                                done
                            done
                        done
                    fi
                done
            done
        done
    done
fi 

if [[ "$1" == "L2" || "$1" == "ALL" ]]; 
then
    #group2: best configurations from group1: change mi_function parameter
    echo "under construction"
fi

if [[ "$1" == "L3" || "$1" == "ALL" ]]; 
then
    #group3: best configuration from group1 applied to larger datasets
    echo "under construction"
fi

if [[ "$1" == "S1" || "$1" == "ALL" ]]; 
then
    #group4: parameters exploration using Gurobi and GLPK
    julia jgvtspsolver.jl gurobi ./input/burma14.xml 10000 secs none
    julia jgvtspsolver.jl gurobi ./input/burma14.xml 10000 combs none
    julia jgvtspsolver.jl gurobi ./input/burma14.xml 10000 both none
fi

if [[ "$1" == "S2" || "$1" == "ALL" ]]; 
then
    #group5: optmization comparison using Gurobi and Christofides
    julia jgvtspsolver.jl gurobi ./input/burma14.xml 10000 both none
    julia jgvtspsolver.jl gurobi ./input/burma14.xml 10000 both christofides,0.00001,2,5pct,lbdecay1,50,0.5,2opt,complementary
    julia jgvtspsolver.jl gurobi ./input/att48.xml 10000 both none
    julia jgvtspsolver.jl gurobi ./input/att48.xml 10000 both christofides,0.00001,2,5pct,lbdecay1,50,0.5,2opt,complementary
    julia jgvtspsolver.jl gurobi ./input/brazil58.xml 10000 both none
    julia jgvtspsolver.jl gurobi ./input/brazil58.xml 10000 both christofides,0.00001,2,5pct,lbdecay1,50,0.5,2opt,complementary
fi

if [[ "$1" == "S3" || "$1" == "ALL" ]]; 
then
    #group6: optmization comparison using Gurobi and Christofides, bigger datasets
    julia jgvtspsolver.jl gurobi ./input/a280.xml 10000 both none
    julia jgvtspsolver.jl gurobi ./input/a280.xml 10000 both christofides,0.00001,2,5pct,lbdecay1,50,0.5,2opt,complementary
    julia jgvtspsolver.jl gurobi ./input/d493.xml 10000 both none
    julia jgvtspsolver.jl gurobi ./input/d493.xml 10000 both christofides,0.00001,2,5pct,lbdecay1,50,0.5,2opt,complementary
fi    