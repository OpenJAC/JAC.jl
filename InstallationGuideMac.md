
## Installation Guide for JAC

The following instructions guide you *stepwise* through the installation of all JAC dependencies.
These instructions have been contributed by Martino Trassinelli (2020) and written for and tested on 
**Mac XX.yy**. In particular, you should be able to run JAC after the following three steps:


**1)** First, we need to install the **julia (version 1.0 or later) programming language**. 
This can be done easily by typing
```
> sudo apt-get install julia
```

- to download the official installer at https://julialang.org/
- create a bin directory in your home if needed: open a terminal and tape 'mkdir bin'


If there is an issue because of an older versions of Ubuntu or some other Linux distributions, 
there is an alternative way to install julia: Download julia installation files from:
**https://julialang.org/downloads/** and extract them to your preferred location.
Afterwards, please, open a terminal and type:
```
> sudo ln -s /FULL/PATH/TO/JULIA/../bin/julia /usr/local/bin/julia
```
in order to create a link to the executable in your system. 
        
        
- if needed put your new bin directory in executable path: add the line 'export PATH=$PATH:$HOME/bin' in .bashrc (or equivalent)
- create a simlink in your $HOME/bin directory: 'ln -s  /Applications/Julia-VERSION.app/Contents/Resources/julia/bin/julia $HOME/bin/.' 
  where VERSION has to be changed with your Julia version you downloaded


**2)** If you don't have a **gfortran compiler** install this by typing
```
> sudo apt-get install gfortran
```
in the terminal. 

Now, you can run julia in a terminal by executing:
```
> julia
```

- for gfortran, use homebrew: 'brew install gcc'  (gfortran is included) or 'brew install gfortran' directly
- enter in Julia  and do what you suggest to do already in the installation file ....

- compile the additional Fortran library for angular coefficients:
- go to the library of interest: from your home directory in a terminal write 
  'cd .julia/packages/JAC/DIRECTORY_NAME_THAT_COULD_CHANGE/deps/ratip2013-angular-coefficients/'  
  (The directory name I think is referring to the name of the last github commit).
- make a local bin file taping  'mkdir ../bin'
- compile the library taping 'make'
- For using JAC as stand alone you do the same in the directory 'JAC/deps/ratip2013-angular-coefficients/'



Press the key "]" to enter the *julia package manager*; you will see the prompt will change 
to pkg>, and type:
```
pkg> add IJulia
```
to **install IJulia**. 

Next, **install the JAC package** by executing:
```
pkg> add https://github.com/OpenJAC/JAC.jl
```
All JAC data files should be now in a folder in "/home/USER/.julia/packages/JAC" 
(please, don't forget to change "USER" in the address by your user name).

Exit julia by typing:
```
> exit()
```


**3)** Install **jupyter environment** by executing the following command in the terminal 
(jupyter can also be installed from the software center)
```
> sudo apt-get install jupyter
```

To run **JAC tutorials**, simply type (don't forget to change "USER" in the address by your user name)
```
> jupyter notebook /home/USER/.julia/packages/JAC
```

This will start a new jupyter session in your default browser. Open a tutorial from the folder 
JAC.jl-master/tutorials.


**Update of the JAC package**

JAC can be updated by running julia, pressing the key "]" and typing
```
pkg> update JAC
```
This will create a new folder in "/home/user/.julia/packages/JAC" with the up-to-date version 
of JAC.



