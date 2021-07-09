
file_data="ddflow.dat";
file_data_tree="ddflow_tree.dat";
file_data_pid="ddflow_pid.dat";

simsets=c("1_default","2_default_robust","3_K_sharp_peak","4_K_flat_top","5_K_bimodal","6_grady_depend_x","7_kernel_platykurtic", "8_kernel_leptokurtic","9_kernel_asymetric","10_mutation_y_rare");

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
    nspe=length(z[[1]]);
    cat(timen,nspe,"\n",file=a$file_data,append=TRUE);
    cat(a$phe$pid+1,"\n",file=a$file_data,append=TRUE);
    for(i in 1:nspe)cat(z$x[i],z$y[i],n[i],"\n",file=a$file_data,append=TRUE);
    cat("\n",file=a$file_data,append=TRUE);
}


edge_y=1000.0;
halt_func <- function(){
    if(max(a$phe$y)>edge_y)a$sparam$flag_halt<<-TRUE;
}


##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##
##_/_/_/_/    Modification of eco-evolutionary settings       _/_/_/_/##
##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##

if(simsets[setid]=="2_default_robust"){
   cat("\nsimulation for ",simsets[setid],"\n");
    x_init=0.0;y_init=0.0;
}

if(simsets[setid]=="3_K_sharp_peak"){
    cat("\nsimulation for ",simsets[setid],"\n");
    x_init=0.0;y_init=0.0;
    po=2;
    sK=10.0;
    K <- function(x){
    return (K0*1.0/(sK*((x)^po)+1.0));
    }
}

if(simsets[setid]=="4_K_flat_top"){
    cat("\nsimulation for ",simsets[setid],"\n");
    x_init=0.0;y_init=0.0;
    po=6;
    sK=1.0;
    K <- function(x){
    return (K0*1.0/(sK*(x^po)+1.0));
    }
}

if(simsets[setid]=="5_K_bimodal"){
    cat("\nsimulation for ",simsets[setid],"\n");
    x_init=0.0;y_init=0.0;
    sK=0.6;sK1=0.6;mmm0=0.0;mmm1=2.0;qqq=0.8;
    K <- function(x){
        return (K0*(exp(-(x-mmm0)^2/(2.0*(sK^2)))+qqq*exp(-(x-mmm1)^2/(2.0*(sK1^2)))));
    }
}


if(simsets[setid]=="6_grady_depend_x"){
    cat("\nsimulation for ",simsets[setid],"\n");
    x_init=0.0;y_init=0.0;
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
}

if(simsets[setid]=="7_kernel_platykurtic"){
    cat("\nsimulation for ",simsets[setid],"\n");
    x_init=0.0;y_init=0.0;
    tail=1.5;
    alpha <- function (x0,x1){
        return (exp(-(x0-x1)^2*(1+tail^2*abs(x0-x1)/sa)/(1+abs(x0-x1)/sa)/(2.0*sa^2)));
    }

    xx=seq(-sa*8,sa*8,,2000);
    sa1=sqrt(sum(alpha(xx,xx*0)*xx^2)/sum(alpha(xx,xx*0)))
    sa=sa^2/sa1;
    sa_adj=sqrt(sum(alpha(xx,xx*0)*xx^2)/sum(alpha(xx,xx*0)));
    cat("tail:", tail," adjusted_sa:",sa," realized_sa:",sa_adj,"\n")
}



if(simsets[setid]=="8_kernel_leptokurtic"){
    cat("\nsimulation for ",simsets[setid],"\n");
    x_init=0.0;y_init=0.0;
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
    x_init=0.0;y_init=0.0;
    bias=0.075;
    alpha <- function (x0,x1){
        return (exp(bias^2/(2.0*sa^2))*exp(-(x0-x1-bias)^2/(2.0*sa^2)));
    }
    
}

if(simsets[setid]=="10_mutation_y_rare"){
    cat("\nsimulation for ",simsets[setid],"\n");
    x_init=0.0;y_init=0.0;
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
}

    
##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##
##_/_/_/_/               Run simulation                       _/_/_/_/##
##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##

n_init=K(x_init);

z=list(x=c(x_init),y=c(y_init));
n=c(n_init);
en=n;

if(0){
arg = commandArgs(trailingOnly=TRUE);
file_data="test.dat";
file_id=NULL;
runid=1;
print(arg);
if(length(arg)>0){
    cat("arguments:\n");
for(i in 1:length(arg)){
    cat(arg[i],"\n");
    eval(parse(text=arg[i]));
}
}
}

runid=1;
set.seed(seed); 
mydir=".simevol/";
