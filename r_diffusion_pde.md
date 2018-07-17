Load Partial Differential Equation Libraries
--------------------------------------------

``` r
suppressWarnings(suppressMessages(library(deSolve)))
suppressWarnings(suppressMessages(require(ReacTran)))
```

Explained 1D Diffusion-Reaction Equation Initial Values and Boundary Conditions
-------------------------------------------------------------------------------

In this example we consider the 1\_D diffusion-reaction model:

$\\frac{\\partial C}{\\partial t}= \\frac{\\partial}{\\partial x} \\left(D \\cdot \\frac{\\partial C}{\\partial x} \\right)-Q$
with *C* the concentration,
*t* the time,
*x* the distance from the origin,
*Q* the consumption rate,
and the following boundary conditions: $\\frac{\\partial C}{\\partial x\_{x=0}} = 0$
*C*<sub>*x* = 10</sub> = *C*<sub>*e**x**t*</sub>
and we create `Grid` as 10 cm (**L**) into 1000 boxes (**N**) e.g. 100mm, or 0.1mm per step

Set up Coefficients & 1D function
---------------------------------

``` r
Grid <-setup.grid.1D(N=1000, L=10)
D <- 1 #diffusion constant
Q <- 1 #uptake rate
Cext <- 20
pde1D <-function(t, C, params){
  tran=tran.1D(C=C, D=D, C.down=Cext, dx= Grid)$dC
  list(tran-Q)
}
```

Model the Diffusion-Reaction Equation Over 100 Time Points (Days)
-----------------------------------------------------------------

``` r
times <- seq(0,100,by=1) #time in days
out <- ode.1D(y=rep(1, Grid$N),times=times,func=pde1D,parms=NULL,nspec=1)
tail(out[,1:11],n=5) # the last five time points, for positions 1:10 ~1mm
```

    ##        time         1         2         3         4         5         6
    ##  [97,]   96 -27.43429 -27.43419 -27.43401 -27.43372 -27.43335 -27.43288
    ##  [98,]   97 -27.49682 -27.49672 -27.49654 -27.49626 -27.49588 -27.49541
    ##  [99,]   98 -27.55783 -27.55773 -27.55754 -27.55726 -27.55689 -27.55642
    ## [100,]   99 -27.61735 -27.61725 -27.61706 -27.61678 -27.61641 -27.61594
    ## [101,]  100 -27.67542 -27.67532 -27.67513 -27.67485 -27.67447 -27.67400
    ##                7         8         9        10
    ##  [97,] -27.43232 -27.43166 -27.43091 -27.43007
    ##  [98,] -27.49485 -27.49419 -27.49344 -27.49260
    ##  [99,] -27.55585 -27.55520 -27.55444 -27.55360
    ## [100,] -27.61537 -27.61471 -27.61396 -27.61311
    ## [101,] -27.67344 -27.67278 -27.67202 -27.67118

Visualize the Diffusion Equation over 100 Days, Approximates Steady-State
-------------------------------------------------------------------------

``` r
image(out, xlab="Time (days)", ylab="Distance (cm)", main="Diffusion PDE", add.contour=TRUE)
```

![](r_diffusion_pde_files/figure-markdown_github/Graph%20the%20100%20Time%20Points-1.png)
