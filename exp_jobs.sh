#!/bin/bash
#experiments to compare performance of christofides and 2factor
julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 none lagrangean
julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 none lagrangean
julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 none lagrangean
julia jgvtspsolver.jl christofides ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 none lagrangean
julia jgvtspsolver.jl factor2approximation ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.5 none lagrangean
julia jgvtspsolver.jl factor2approximation ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.5 none lagrangean
julia jgvtspsolver.jl factor2approximation ./input/burma14.xml 10000 0.00001 2 1pct lbdecay1 50 0.75 none lagrangean
julia jgvtspsolver.jl factor2approximation ./input/burma14.xml 10000 0.00001 2 1pct lbdecay2 50 0.75 none lagrangean
