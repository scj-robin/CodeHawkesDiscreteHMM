This repo contains an R implementation of the code associated with the paper "Markov switching discrete-time Hawkes process: application to the monitoring of bats' behavior", A. Bonnet, S. Robin (2025). The data is extracted from the publicly availabe dataset provided by the Vigie-Chiro program: https://www.vigienature.fr/fr/chauves-souris




### Data

$n = 524$ calls along one night (times have been rescaled to the (0, 1) interval).


```{r data, echo=FALSE}
setwd('./')
times <- read.table('batcalls.csv', header=FALSE)[, 1]
n <- length(times)
plot(times, 1:n, type='b', pch=20, xlab='time', ylab='count')
cat('Number of events: n =', n, '\n')
```
<img src="./plot_trajectory.png" width="500">

```


### Libraries and fixing parameters

```{r dims, echo=FALSE}
library('hawkesbow') 
source('Functions.R')
Tmax <- 1
Q <- 3 # number of hidden states
coef <- 2 # discretisation coefficient: N = coef * n
N <- ceiling(coef*n)
cat('Tmax =', Tmax, ', coef =', coef, ', N =', N, '\n')



### Builds the discrete-time path

```{r discrete, echo=FALSE}
contPath <- list(times=times, states=rep(1, n), hiddenPath=list(states=1, times=0))
discPath <- MakeDiscreteHawkes1D(contPath=contPath, N=N, Tmax=Tmax)
cat('Observed counts : \n', 
    discPath$y[1:20], '\n')
```

### Initialization

```{r init, echo=FALSE}
initFile <- paste0('batcalls-init-Q', Q, '.Rdata')
if(!file.exists(initFile)){
  init <- InitEMHMMHawkes1D(contPath=contPath, nInit=N, Q=Q, Tmax=Tmax)
  save(init, file=initFile)
}else{load(initFile)}
```

### Fit the discrete-time Hawkes HMM

```{r fit, echo=FALSE}
fitFile <- paste0('batcalls-fit-Q', Q, '.Rdata')
if(!file.exists(fitFile)){
  fit <- EMHMMHawkes1D(discPath=discPath, discParms=init, Q=Q)
  save(fit, file=fitFile)
}else{load(fitFile)}
cat('Parms : \n')
cat('alpha =', fit$discParms$alpha, ', beta =', fit$discParms$beta, '\n', 
    'mu =', fit$discParms$mu, '\n', 
    'diag(pi) =', diag(fit$discParms$pi), '\n', 
    'state frequencies =', table(fit$classif)/N, '\n', 
    'logLik = ', fit$logL, '\n')
```

### Output

```{r plot, echo=FALSE}
plot((1:N)/N, cumsum(discPath$y), type='b', pch=20, col=fit$classif, xlab='time', ylab='count')
abline(v=(0.5+which(diff(fit$classif)!=0))/N, lty=2, lwd=2)
```
<img src="./segmented_trajectory.png" width="500">

