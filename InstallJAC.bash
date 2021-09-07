#!/bin/bash
tput bold && tput setaf 6; echo 'Downloading and installing julia'
tput sgr0   
sudo wget https://julialang-s3.julialang.org/bin/linux/x64/1.6/julia-1.6.2-linux-x86_64.tar.gz
sudo tar -xvzf julia-1.6.2-linux-x86_64.tar.gz -C /usr/local
sudo ln -s /usr/local/julia-1.6.2/bin /usr/local/bin/julia
sudo rm julia-1.6.2-linux-x86_64.tar.gz
tput bold && tput setaf 6; echo 'Installing jupyter'
tput sgr0   
sudo apt-get install jupyter
tput bold && tput setaf 6; echo 'Installing JAC and its dependencies'
tput sgr0   
julia JuliaSetup.jl
