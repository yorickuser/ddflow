## file ddflow_osm0d.R
## R program for simulation of DD-flows with OSM
## (oligomorphic stochastic model, Ito and Dieckmann (2007, 2014)),
## written by Hiroshi C. Ito 2021.
## OSM was implemented following Ito and Sasaki (2020) as an R-package "simevol" (no documentation yet).

## simevol can be downloaded (from my github account "yorickuser") and installed by using devtools.
## library(devtools)
## install_github('yorickuser/simevol@0.1.2')
##
## Eco-evolutionary setting can be chosen from the list in "simsets" (e.g., "3_K_sharp_peak" is selected by "setid=3").
##
## Execution method 1 (at R prompt): source("ddflow_osm0d.R")
## Execution method 2 (at console): Rscript ddflow_osm0d.R setid=3



## copyright (C) 2021 Hiroshi C. Ito
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License Version 2 as 
## published by the Free Software Foundation.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## A copy of the GNU General Public License is available at
## http://www.r-project.org/Licenses/


library(simevol);

ofile="ddflow";

simsets=c("1_default",
          "2_default_robust",
          "3_K_sharp_peak",
          "4_K_flat_top",
          "5_K_bimodal",
          "6_grady_depend_x",
          "7_kernel_platykurtic",
          "8_kernel_leptokurtic",
          "9_kernel_asymetric",
          "10_mutation_y_rare",
          "11_multiple_geographic_regions");

runnames=c("1_def","2_defr","3_Ksharp","4_Kflat","5_Kbimod","6_gradyx","7_plat","8_lept","9_asym","10_yrare","11_geog");

setid=1;

##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##
##_/_/_/_/    Default eco-evolutionary settings           _/_/_/_/##
##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##

K0=1.0; ##maximum carrying capacity
sK=1.0; ##sd of carrying capacity
sa=0.15; ##sd of competition kernel
sm=0.01; ##sd of mutation
grady=0.1; ## fitness gradient for fundamental trait
m_rate=1.0; ## mutation rate per unit population size
seed=13; ## rundam seed
x_init=-0.2; 
y_init=0.0;


K <- function(x){
    return (K0*exp(-(x^2)/(2.0*sK^2)));
}


mutate <- function(z){
        mx=z$x+rnorm(1,mean=0.0,sd=sm);
        my=z$y+rnorm(1,mean=0.0,sd=sm);    
     return(list(x=mx,y=my));
}

alpha <- function (x0,x1){
    return (exp(-((x0-x1)^2)/(2.0*sa^2)));
}


fitness <-function(z1,z,n){
    return(grady*(z1$y)+1.0-sum(alpha(z$x,z1$x)*n)/K(z1$x));
}

set_parms <-function(z,n){
    Kp=K(z$x);
    ap=matrix(n,nrow=length(n),ncol=length(n),byrow=T)*0;
    for(i in 1:length(n)){
        for(j in 1:length(n)){
            ap[i,j]=alpha(z$x[j],z$x[i]);
        }
    }
    return(list(z=z,Kp=Kp,ap=ap));
}

pop_dynamics <- function(t,n,parms){
    list(n*(grady*parms$z$y+1-rowSums(parms$ap%*%n)/parms$Kp));     
}


output <- function(timen,z,n){
    options(scipen=100);
    nspe=length(z[[1]]);
    cat(timen,nspe,"\n",file=a$file_data,append=TRUE);
    cat(a$phe$pid+1,"\n",file=a$file_data,append=TRUE);
    for(i in 1:nspe)cat(z$x[i],z$y[i],n[i],"\n",file=a$file_data,append=TRUE);
    cat("\n",file=a$file_data,append=TRUE);
    options(scipen=0);
}


edge_y=1000.0;
halt_func <- function(){
    if(max(a$phe$y)>edge_y)a$sparam$flag_halt<<-TRUE;
}

z=list(x=c(x_init),y=c(y_init));
n_init=K(x_init);
n=c(n_init);en=n;

##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##
##_/_/_/_/    Modification of eco-evolutionary settings       _/_/_/_/##
##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##


arg = commandArgs(trailingOnly=TRUE);
runid=1;
show_comwin=0;
bgid=2;
print(arg);
if(length(arg)>0){
    cat("arguments:\n");
for(i in 1:length(arg)){
    cat(arg[i],"\n");
    eval(parse(text=arg[i]));
}
}



if(simsets[setid]=="2_default_robust"){
   cat("\nsimulation for ",simsets[setid],"\n");
   x_init=0.0; y_init=0.0; z=list(x=c(x_init),y=c(y_init));
   n_init=K(x_init);
   n=c(n_init);en=n;
}

if(simsets[setid]=="3_K_sharp_peak"){
    cat("\nsimulation for ",simsets[setid],"\n");
    x_init=0.0;y_init=0.0; z=list(x=c(x_init),y=c(y_init));
    po=2;
    sK=10.0;
    K <- function(x){
    return (K0*1.0/(sK*((x)^po)+1.0));
    }
    n_init=K(x_init);
    n=c(n_init);en=n;
}

if(simsets[setid]=="4_K_flat_top"){
    cat("\nsimulation for ",simsets[setid],"\n");
    x_init=0.0;y_init=0.0; z=list(x=c(x_init),y=c(y_init));
    po=6;
    sK=1.0;
    K <- function(x){
    return (K0*1.0/(sK*(x^po)+1.0));
    }
    n_init=K(x_init);
    n=c(n_init);en=n;
}

if(simsets[setid]=="5_K_bimodal"){
    cat("\nsimulation for ",simsets[setid],"\n");
    x_init=0.0;y_init=0.0; z=list(x=c(x_init),y=c(y_init));
    sK=0.6;sK1=0.6;mmm0=0.0;mmm1=2.0;qqq=0.8;
    K <- function(x){
        return (K0*(exp(-(x-mmm0)^2/(2.0*(sK^2)))+qqq*exp(-(x-mmm1)^2/(2.0*(sK1^2)))));
    }
    n_init=K(x_init);
    n=c(n_init);en=n;
}


if(simsets[setid]=="6_grady_depend_x"){
    cat("\nsimulation for ",simsets[setid],"\n");
    x_init=0.0;y_init=0.0; z=list(x=c(x_init),y=c(y_init));
    offset_grady=1.0;
    sd_grady=0.55;
    m_grady=1.5;
    gradyy <-function(x){
        return (grady*(offset_grady+ exp(-(x-m_grady)*(x-m_grady)/(2*sd_grady*sd_grady))));
    }
    
    fitness <-function(z1,z,n){
        return(gradyy(z1$x)*(z1$y)+1.0-sum(alpha(z$x,z1$x)*n)/K(z1$x));
    }
    
    pop_dynamics <- function(t,n,parms){
        list(n*(gradyy(parms$z$x)*parms$z$y+1-rowSums(parms$ap%*%n)/parms$Kp));     
    }
    n_init=K(x_init);
    n=c(n_init);en=n;

}

if(simsets[setid]=="7_kernel_platykurtic"){
    cat("\nsimulation for ",simsets[setid],"\n");
    x_init=0.0;y_init=0.0; z=list(x=c(x_init),y=c(y_init));
    tail=1.5;
    alpha <- function (x0,x1){
        return (exp(-(x0-x1)^2*(1+tail^2*abs(x0-x1)/sa)/(1+abs(x0-x1)/sa)/(2.0*sa^2)));
    }

    xx=seq(-sa*8,sa*8,,2000);
    sa1=sqrt(sum(alpha(xx,xx*0)*xx^2)/sum(alpha(xx,xx*0)))
    sa=sa^2/sa1;
    sa_adj=sqrt(sum(alpha(xx,xx*0)*xx^2)/sum(alpha(xx,xx*0)));
    cat("tail:", tail," adjusted_sa:",sa," realized_sa:",sa_adj,"\n")

    n_init=K(x_init);
    n=c(n_init);en=n;

}



if(simsets[setid]=="8_kernel_leptokurtic"){
    cat("\nsimulation for ",simsets[setid],"\n");
    x_init=0.0;y_init=0.0; z=list(x=c(x_init),y=c(y_init));
    tail=1.0/1.5;
    alpha <- function (x0,x1){
            return (exp(-(x0-x1)^2*(1+tail^2*abs(x0-x1)/sa)/(1+abs(x0-x1)/sa)/(2.0*sa^2)));
    }

    xx=seq(-sa*8,sa*8,,2000);
    sa1=sqrt(sum(alpha(xx,xx*0)*xx^2)/sum(alpha(xx,xx*0)))
    sa=sa^2/sa1;
    sa_adj=sqrt(sum(alpha(xx,xx*0)*xx^2)/sum(alpha(xx,xx*0)));
    cat("tail:", tail," adjusted_sa:",sa," realized_sa:",sa_adj,"\n")

}

if(simsets[setid]=="9_kernel_asymetric"){
    cat("\nsimulation for ",simsets[setid],"\n");
    x_init=0.0;y_init=0.0; z=list(x=c(x_init),y=c(y_init));
    bias=0.075;
    alpha <- function (x0,x1){
        return (exp(bias^2/(2.0*sa^2))*exp(-(x0-x1-bias)^2/(2.0*sa^2)));
    }
    n_init=K(x_init);
    n=c(n_init);en=n;

}

if(simsets[setid]=="10_mutation_y_rare"){
    cat("\nsimulation for ",simsets[setid],"\n");
    x_init=0.0;y_init=0.0; z=list(x=c(x_init),y=c(y_init));
    grady=0.025;
    m_rate_x=1.0;
    sm_x=0.01;
    sm_y=2.0;
    m_rate_y=m_rate_x*sm_x^2/sm_y^2;
    m_rate=m_rate_x+m_rate_y;

    mutate <- function(z){
        if(runif(1)<m_rate_x/(m_rate_x+m_rate_y)){
            mx=z$x+rnorm(1,mean=0.0,sd=sm_x);
            my=z$y;
        }else{
            mx=z$x;
            my=z$y+rnorm(1,mean=0.0,sd=sm_y);    
        }
        return(list(x=mx,y=my));
    }
    n_init=K(x_init);
    n=c(n_init);en=n;

}

if(simsets[setid]=="11_multiple_geographic_regions"){
    cat("\nsimulation for ",simsets[setid],"\n");
    x_init=0.0;y_init=0.0;w_init=0.0;
    z=list(x=c(x_init),y=c(y_init),w=c(w_init));


    sKg=0.1;
    sag=0.0001;
    mrate_trait=1.0;
    mrate_geog=0.00005;
    m_rate=mrate_trait+mrate_geog;
    smg=0.1;
    seed=17;
    bgid=3;
    K <- function(x,w){
        return (K0*exp(-(x^2)/(2.0*sK^2))*exp(-(w^2)/(2.0*sKg^2)));
    }

    n_init=K(x_init,w_init);
    n=c(n_init);en=n;

    mutate <- function(z){
        if(runif(1)<mrate_geog/m_rate){
            if(runif(1)<0.5){
                mw=z$w+smg;
            }
            else{
                mw=z$w-smg;
            }
            mx=z$x;
            my=z$y;
        }
        else{
            mx=z$x+rnorm(1,mean=0.0,sd=sm);
            my=z$y+rnorm(1,mean=0.0,sd=sm);
            mw=z$w;
        }
        return(list(x=mx,y=my,w=mw));
    }
    
    alpha <- function (x0,w0,x1,w1){
        return (exp(-((x0-x1)^2)/(2.0*sa^2))*as.numeric(abs(w0-w1)<0.0001));
    }


    fitness <-function(z1,z,n){
        return(grady*(z1$y)+1.0-sum(alpha(z$x,z$w,z1$x,z1$w)*n)/K(z1$x,z1$w));
    }
    
    set_parms <-function(z,n){
        Kp=K(z$x,z$w);
        ap=matrix(n,nrow=length(n),ncol=length(n),byrow=T)*0;
        for(i in 1:length(n)){
            for(j in 1:length(n)){
                ap[i,j]=alpha(z$x[j],z$w[j],z$x[i],z$w[i]);
            }
        }
        return(list(z=z,Kp=Kp,ap=ap));
    }
    
    pop_dynamics <- function(t,n,parms){
        list(n*(grady*parms$z$y+1-rowSums(parms$ap%*%n)/parms$Kp));     
    }
    
    output <- function(timen,z,n){
        options(scipen=100);
        nspe=length(z[[1]]);
        cat(timen,nspe,2,"\n",file=a$file_data,append=TRUE);
        cat(a$phe$pid+1,"\n",file=a$file_data,append=TRUE);
        for(i in 1:nspe)cat(z$x[i],z$w[i],z$y[i],n[i],"\n",file=a$file_data,append=TRUE);
        cat("\n",file=a$file_data,append=TRUE);
        options(scipen=0);        
    }
    
}

##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##
##_/_/_/_/               Run simulation                       _/_/_/_/##
##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##

ofile=paste0(ofile,"_",runnames[setid]);
file_data=paste0(ofile,".dat");
file_data_tree=paste0(ofile,"_tree.dat");
file_data_pid=paste0(ofile,"_pid.dat");
cat("\n output files:\n");
print(file_data);
print(file_data_tree);
print(file_data_pid);

set.seed(seed); 
mydir=".simevol/";
if(show_comwin==1)comwin(1,mydir=mydir);

##.ee.set();
simevol(z,en,fitness=fitness,
        mutate=mutate,
        fitness_contour=FALSE,
        show_interval=80,
        out_interval=80,
        pop_dynamics=pop_dynamics,
        set_parms=set_parms,
        output=output,
        halt_func=halt_func,
        file_data=file_data,
        file_data_tree=file_data_tree,
        file_data_pid=file_data_pid,
        m_rate=m_rate,
        runid=runid,
        runname=runnames[setid],
        show_subwin=TRUE,
        amp_invf_fix=TRUE,
        amp_invf=1.0,
        mydir=mydir,
        bgid=bgid,
        edge_extinct=1e-6,
        n_mutant_init=1e-5);
