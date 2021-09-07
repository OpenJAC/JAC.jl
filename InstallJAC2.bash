#!/bin/bash
tput bold && tput setaf 6; echo 'Downloading and installing julia'
tput sgr0   
sudo snap install julia --classic
tput bold && tput setaf 6; echo 'Installing jupyter'
tput sgr0   
sudo apt-get install jupyter
tput bold && tput setaf 6; echo 'Installing JAC and its dependencies'
tput sgr0   
julia JuliaSetup.jl
