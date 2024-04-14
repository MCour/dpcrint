---
title: "Using-dpcrint"
output: 
  html_document:
    keep_md: true
---




```r
library(dpcrint)
par(mar=c(3,3,1,1))
```

## Genome integrity estimation from multiplex dPCR

This package was built in order to recover, from multiplex dPCR, an idea of the integrity of the genome in a sample. In this kind of study, $J$ positions along the genome are tageted by a multiplex dPCR.
The presence or absence of signal for these targets gives an idea of the integrity of the genome in the sample analysed. 

A first approach to carry out this analysis is to look directly at the percentage of partitions showing a signal for all the targets (among the partitions for which a signal is observed). 
However, this does not allow the following possibilities to be taken into account: 

* among the partitions showing signals for all $J$ targets, the presence of several incomplete sequences, each covering part of the targets,
* among the partitions showing signals for only part of the targets, the presence of complete sequences for which PCR amplification did not work for one or more of the targets. 

In addition, the integrity of the sequences may be partial and it would be useful to evaluate it too.

These reasons led us to develop a deconvolution algorithm adapted to these data. It takes as inputs: 

* For each combination of signals present/absent for each target, the number of concerned partitions,
* the expected average number of sequences per partition (supplied by the manufacturer of the dPCR device),
* the probability of amplification for each target on the genome. 

If some part of this information is missing, it can be recovered at some conditions.

We here illustrate the implemented methods on simulated example data sets.


## Recovering the amplification probabilities

If the probability of amplification for each target on the genome is not known, a preliminary rapid estimate of those probabilities can be made based on experiments carried out with samples containing only complete sequences.

To simulate such data set, one can use the **SimulObs** function:

```r
N=26000 # number of partitions
lambda=0.05 # average number of sequences in one partition
J=5 # number of targets
pobs=c(0.8,0.6,0.75,0.99,0.98) # amplification probability for each target
coverseq=defseqtypes(J) # sequence types definition in terms of target coverage
L=nrow(coverseq)
propl=numeric(L) # proportion of each sequence type
propl[1]=1 
# print(cbind(coverseq,propl)) 
ressim=SimulObs(N,lambda,propl,pobs)
ObsTypes=ressim$ObsTypes
NY=ressim$NY # number of each observation type in terms of target coverage
print(NY)
#> Yobs
#> 0-0-0-0-0 0-0-0-1-0 0-0-0-1-1 0-0-1-1-0 0-0-1-1-1 0-1-0-0-1 0-1-0-1-0 0-1-0-1-1 
#>     24719         2        20         1        74         2         2        34 
#> 0-1-1-0-1 0-1-1-1-1 1-0-0-0-1 1-0-0-1-0 1-0-0-1-1 1-0-1-0-1 1-0-1-1-0 1-0-1-1-1 
#>         1       112         1         2        87         4         9       304 
#> 1-1-0-0-1 1-1-0-1-0 1-1-0-1-1 1-1-1-0-1 1-1-1-1-0 1-1-1-1-1 
#>         1         3       155         5         6       456
# identical(coverseq,ressim$coverseq)
```

To recover the amplification probability of each target from such data, one can use the **DirectOptPobs** function:

```r
respobs=DirectOptPobs(NY,ObsTypes,lambda0=lambda)
print(respobs)
#> [1] 0.8025311 0.6006433 0.7541587 0.9887777 0.9800028
print(pobs)
#> [1] 0.80 0.60 0.75 0.99 0.98
```

To recoverthe amplification probability of each target from such data together with the average number of sequences in each partition, one can use the **DirectOptPobs** function with the option *evallam=TRUE* :

```r
respobslam=DirectOptPobs(NY,ObsTypes,evallam=TRUE)
print(respobslam$pobs)
#> [1] 0.8020107 0.5982952 0.7531460 0.9881764 0.9796605
print(pobs)
#> [1] 0.80 0.60 0.75 0.99 0.98
print(respobslam$lambda)
#> [1] 0.05052566
print(lambda)
#> [1] 0.05
```

<!-- One can also try to recover those parameters from other dPCR results from samples with known sequence types proportions through the SEM algorithm, but the results that are obtained can be very poor. -->

<!-- ```{r datasim0} -->
<!-- respobsSEM=OptPobsSEM(NY,ObsTypes,propl,coverseq,lambda0=lambda,evallam=F,Nit=100,plotopt=T) -->
<!-- print(respobsSEM$pobs) -->
<!-- print(pobs) -->
<!-- ``` -->

## Recovering the sequence types proportions in sample

In order to recover the sequence types proportions in the sample (from any sample type), one can use the proposed EM algorithm, embedded in the **OptProplEM** function:

```r
respropl=OptProplEM(NY,ObsTypes,pobs=pobs,lambda0=lambda)
abline(h=propl,col=1:L,lty=2)
```

![](/private/var/folders/4h/b18j_vj51rn90w8rzq1hr1bh0000gp/T/RtmpnorXiq/preview-af847bad1679.dir/Using-dpcrint_files/figure-html/optpropl0-1.png)<!-- -->

```r
plot(as.vector(respropl$propl)~propl,ylab="estimated proportions",xlab="true proportions",pch=19)
abline(0,1,lty=2)
```

![](/private/var/folders/4h/b18j_vj51rn90w8rzq1hr1bh0000gp/T/RtmpnorXiq/preview-af847bad1679.dir/Using-dpcrint_files/figure-html/optpropl0-2.png)<!-- -->

In order to recover both the sequence types proportions in the sample (from any sample type) and the average number of sequences in one partition, one can use the proposed EM algorithm with the *evallam=TRUE* option:

```r
respropllam=OptProplEM(NY,ObsTypes,pobs=pobs,evallam=T)
abline(h=propl,col=1:L,lty=2)
```

![](/private/var/folders/4h/b18j_vj51rn90w8rzq1hr1bh0000gp/T/RtmpnorXiq/preview-af847bad1679.dir/Using-dpcrint_files/figure-html/optpropllam0-1.png)<!-- -->

```r
plot(as.vector(respropllam$propl)~propl,ylab="estimated proportions",xlab="true proportions",pch=19)
abline(0,1,lty=2)
```

![](/private/var/folders/4h/b18j_vj51rn90w8rzq1hr1bh0000gp/T/RtmpnorXiq/preview-af847bad1679.dir/Using-dpcrint_files/figure-html/optpropllam0-2.png)<!-- -->

```r
plot(respropllam$LAMBDA,ylab="lambda",xlab="iterations")
abline(h=lambda)
```

![](/private/var/folders/4h/b18j_vj51rn90w8rzq1hr1bh0000gp/T/RtmpnorXiq/preview-af847bad1679.dir/Using-dpcrint_files/figure-html/optpropllam0-3.png)<!-- -->


With another toy data set:

```r
propl=exp(0.99)
for(j in 2:J){
  propl=c(propl,rep(exp(0.99*j),j))
}
propl=propl/sum(propl) # proportion of each sequence type
# cbind(coverseq,propl) 
ressim=SimulObs(N,lambda,propl,pobs)
ObsTypes=ressim$ObsTypes
NY=ressim$NY # number of each observation type in terms of target coverage
print(NY)
#> Yobs
#> 0-0-0-0-0 0-0-0-0-1 0-0-0-1-0 0-0-0-1-1 0-0-1-0-0 0-0-1-1-0 0-0-1-1-1 0-1-0-0-0 
#>     24859       167       210        82       164        60        19       132 
#> 0-1-0-0-1 0-1-0-1-0 0-1-0-1-1 0-1-1-0-0 0-1-1-1-0 0-1-1-1-1 1-0-0-0-0 1-0-0-0-1 
#>         2         3         3        33        14         5       170         2 
#> 1-0-1-0-0 1-0-1-1-0 1-0-1-1-1 1-1-0-0-0 1-1-0-1-0 1-1-1-0-0 1-1-1-1-0 1-1-1-1-1 
#>        13         3         1        39         2        11         5         1
# identical(coverseq,ressim$coverseq)
```


```r
respropl=OptProplEM(NY,ObsTypes,pobs=pobs,lambda0=lambda)
abline(h=propl,col=1:L,lty=2)
```

![](/private/var/folders/4h/b18j_vj51rn90w8rzq1hr1bh0000gp/T/RtmpnorXiq/preview-af847bad1679.dir/Using-dpcrint_files/figure-html/optpropl-1.png)<!-- -->

```r
plot(as.vector(respropl$propl)~propl,ylab="estimated proportions",xlab="true proportions",pch=19)
abline(0,1,lty=2)
```

![](/private/var/folders/4h/b18j_vj51rn90w8rzq1hr1bh0000gp/T/RtmpnorXiq/preview-af847bad1679.dir/Using-dpcrint_files/figure-html/optpropl-2.png)<!-- -->


```r
respropllam=OptProplEM(NY,ObsTypes,pobs=pobs,evallam=T)
abline(h=propl,col=1:L,lty=2)
```

![](/private/var/folders/4h/b18j_vj51rn90w8rzq1hr1bh0000gp/T/RtmpnorXiq/preview-af847bad1679.dir/Using-dpcrint_files/figure-html/optpropllam-1.png)<!-- -->

```r
plot(as.vector(respropllam$propl)~propl,ylab="estimated proportions",xlab="true proportions",pch=19)
abline(0,1,lty=2)
```

![](/private/var/folders/4h/b18j_vj51rn90w8rzq1hr1bh0000gp/T/RtmpnorXiq/preview-af847bad1679.dir/Using-dpcrint_files/figure-html/optpropllam-2.png)<!-- -->

```r
plot(respropllam$LAMBDA,ylab="lambda",xlab="iterations")
abline(h=lambda)
```

![](/private/var/folders/4h/b18j_vj51rn90w8rzq1hr1bh0000gp/T/RtmpnorXiq/preview-af847bad1679.dir/Using-dpcrint_files/figure-html/optpropllam-3.png)<!-- -->
