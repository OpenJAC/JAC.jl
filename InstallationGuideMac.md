
## Installation Guide for JAC

The following instructions guide you *stepwise* through the installation of all JAC dependencies.
These instructions have been contributed by Martino Trassinelli (2020) and written for and tested on
**Mac XX.yy**. In particular, you should be able to run JAC after the following three steps:


**1)** First, we need to install the **julia (version 1.0 or later) programming language**.
This can be done easily
- by typing
```
> sudo apt-get install julia
```
if you use apt-get,
- by typing
```
> brew install julia
```
if you use homebrew,
- by downloading the official installer at https://julialang.org/
(this is the only option fully tested so far)

If you install Julia from the official download, you need to make available the executable in terminal. To do that, proceed with the following steps.  
- if needed create the bin directory in your home. For this open a terminal and tape 'mkdir bin'    
- if needed put your new bin directory in executable path: add the line 'export PATH=$PATH:$HOME/bin' in .bashrc (or equivalent)
- create a simlink in your $HOME/bin directory: 'ln -s  /Applications/Julia-VERSION.app/Contents/Resources/julia/bin/julia $HOME/bin/.'
  where VERSION has to be changed with your Julia version you downloaded

Now, you can run julia in a terminal by executing:
```
> julia
```

**2)** If you don't have a **gfortran compiler** install this
- by typing
```
> sudo apt-get install gfortran
```
in the terminal, if you use apt-get,
- by typing
```
> brew install gfortran
```
or
```
> brew install gcc
```
in the terminal, if you use homebrew (gfortran is inclued in the gcc installation).

To test if gfortran is working, test it in a terminal
```
> gfortran --version
```
You should obtain an output like
```
GNU Fortran (Homebrew GCC 9.2.0_3) 9.2.0
Copyright (C) 2019 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
```


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

**2bis)** Once you downloaded JAC.jl (as package via Julia, or github or direct dowload), you have to compile the additional Fortran library for angular coefficients.

To do it
- go to the library of interest: from your home directory in a terminal write

  ```
  > cd .julia/packages/JAC/DIRECTORY_NAME_THAT_COULD_CHANGE/deps/ratip2013-angular-coefficients
  ```
  The directory name I think is referring to the name of the last github commit.
- make a local bin file taping
```
> mkdir ../bin
```  
'mkdir ../bin'
-  compile the library taping
```
> make
```
- For using JAC as stand alone you do the same in the directory 'JAC/deps/ratip2013-angular-coefficients/'


**3)** Install **jupyter environment** by executing the following command in the terminal
(jupyter can also be installed from the software center)
```
> sudo apt-get install jupyter
```
or
```
> brew install jupyter
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
