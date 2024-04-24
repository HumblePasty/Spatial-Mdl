---
  title: "Gaussian Processes models for point-referenced spatial data, part 1"
subtitle: "Spatial Statistics - BIOSTAT 696/896"
author: "Michele Peruzzi"
institute: "University of Michigan"
execute:
  cache: true
format:
  revealjs:
  auto-stretch: true
center: true
html-math-method: katex
transition: none
slide-number: true
callout-appearance: minimal
code-line-numbers: true
fontsize: 16pt
output-file: l04_gaussianprocesses.html
theme: [simple, custom.css] #, web/pulse_colors.css] # alternative themes (subset): default, night, dark
embed-resources: false
echo: true
fig-dpi: 150
# incremental: true  # Remove comment if you like incremental bullet points
footer: "[Home page](index.html)"
knitr:
  opts_chunk:
  out.width: "60%"
editor: 
  markdown: 
  wrap: 72
editor_options: 
  chunk_output_type: console
---
  
  ## Introducing Gaussian Processes
  
  -   We are considering the case of point-referenced data
-   We want to model spatial dependence of a set of random variables
-   We have one random variable for each location in
$\cal L = \{s_1, \dots, s_n \}$:
  $$\{ Y(s_1), Y(s_2), \dots, Y(s_n) \} = \{ Y(s) : s \in \cal L \}$$
  -   We need to define $p(Y(s_1), \dots, Y(s_n))$, that is our model of
dependence for this specific set of random variables
-   [What about $p(Y(s_1), \dots, Y(s_n)) = \prod_{i=1}^n p(Y(s_i))$
       ?]{.orange}

## Introducing Gaussian Processes

-   Independence assumption (IND):
  $p(Y(s_1), \dots, Y(s_n)) = \prod_{i=1}^n p(Y(s_i))$
  -   IND is very common when we have *replicated* data, i.e. multiple
observations of the same random variable
-   IND is inappropriate with spatial data because it assumes there is
no dependence at all! We want dependence between random variables
that is based on the notion of proximity or distance
-   *The closer in space* $Y(s)$ and $Y(s')$ are, the more they should
    be dependent on each other
-   To model spatial data is to model (spatial) dependence within a set
    of random variables
-   Multivariate Gaussian assumption (MVN): assume
    $p(Y(s_1), \dots, Y(s_n)) = MVN(Y_{\cal L}; \mu, \Sigma)$
-   In MVN, $Y_{\cal L}$ is the column vector
    $Y_{\cal L} = (Y(s_1), \dots, Y(s_n))^\top$
-   Let's write $Y = Y_{\cal L}$ in short to keep notation light
                                        -   $\mu$ and $\Sigma$ are the mean vector and covariance matrix for the
                                        random *vector* $Y$
                                          -   [What is $Cov(Y(s_5), Y(s_8))$ under this model?]{.orange}
                                        
                                        ## Introducing Gaussian Processes
                                        
                                        -   [What is $Cov(Y(s_5), Y(s_8))$ under this model?]{.orange}
                                        $\Sigma_{5,8}$
                                          -   We usually assume $E(Y) = \mu = 0$
                                          -   We have just defined $Cov(Y) = \Sigma$
                                          -   Under our Gaussian assumption, [what else do we need to define to
                                                                              complete our model for $p(Y)$?]{.orange}
                                        
                                        ## Introducing Gaussian Processes
                                        
                                        -   Under our Gaussian assumption, [what else do we need to define to
                                                                            complete our model for $p(Y)$?]{.orange}
                                        -   Nothing. A MVN is completely defined by its mean vector and
                                        covariance matrix
                                        -   [How many parameters does this model have?]{.orange}
                                        
                                        $$ Y \sim MVN(\mu, \Sigma) $$
                                          
                                          ## Introducing Gaussian Processes
                                          
                                          -   Let $\mu = 0$ and estimating $\Sigma$ we have $n(n+1)/2$ parameters
                                        in this model
                                        
                                        $$ Y \sim MVN(0, \Sigma) $$
                                          
                                          -   We estimate $\Sigma_{ij} = Cov(Y(s_i), Y(s_j))$ for all
                                        $i=1,\dots,n$, $j=1,\dots,i$ since $\Sigma$ is symmetric
                                        -   This is a **difficult** problem: very high dimensional and we need
                                        $\Sigma$ positive definite!
                                          
                                          In this model:
                                          
                                          -   The set of locations $\cal L$ is given and fixed
                                        -   Our estimate for $\Sigma$ is dependent on $\cal L$ being exactly
                                        what it is
                                        
                                        [Therefore, is the above model *useful*?]{.orange}
                                        
                                        -   [What if we had a *different* fixed set of locations? How would we
                                             model that?]{.orange}
                                        -   [In other words, if we estimated the model above on $\cal L$, can we
                                             say anything about]{.orange}
                                        
                                        $$ Cov(Y(s_{n+1}), Y(s_{n+2})) $$ [for a pair of locations
                                                                           $s_{n+1}, s_{n+2} \notin \cal L$?]{.orange}
                                        
                                        ## Introducing Gaussian Processes
                                        
                                        -   [What if we had a *different* fixed set of locations? How would we
                                             model that?]{.orange}
                                        -   [In other words, if we estimated the model above on $\cal L$, can we
                                             say anything about]{.orange}
                                        
                                        $$ Cov(Y(s_{n+1}), Y(s_{n+2})) $$ [for a pair of locations
                                                                           $s_{n+1}, s_{n+2} \notin \cal L$?]{.orange}
                                        
                                        -   Even if we were able to estimate that difficult-to-estimate model,
                                        we still would not be able to say anything about
                                        $Cov(Y(s_{n+1}), Y(s_{n+2}))$
                                          
                                          ## Introducing Gaussian Processes
                                          
                                          -   We need a way to learn about $Y(\cdot)$ at any point in the domain,
                                        based on observing data on any finite set $\cal L$
                                          
                                          **Building a Gaussian Process** (overview)
                                        
                                        1)  For any fixed $\cal L$, we assume that $Y_{\cal L}$ is
**multivariate Gaussian**
  2)  We define the finite dimensional distribution for any set of random
variables $Y(s_1), \dots, Y(s_n)$ based on a choice of **covariance
function**
  3)  (Check Kolmogorov's conditions for the existence of a corresponding
    stochastic process: that stochastic process is the Gaussian Process)

## Building Gaussian Processes

-   Consider a set of locations in our spatial domain
    $\cal L \subset \cal D \subset \Re^d$
-   Define the mean function $\mu(\cdot) : \cal D \to \Re$
-   Typical assumption: $\mu(s) = 0\quad \forall s\in \cal D$
-   Define the covariance function
    $C(\cdot, \cdot): \cal D \times \cal D \to \Re$ as a symmetric and positive
    definite function
-   We then let

$$ Y_{\cal L} \sim N(0, C_{\cal L})  $$

-   Where the $(i,j)$ element of the matrix $C_{\cal L}$ is
    $C(s_i, s_j)$
-   [Therefore, $E(Y(s_i) Y(s_j)) =$ ??]{.orange}

## Building Gaussian Processes

-   By construction we have $E(Y(s_i))=0$ and $E(Y(s_j)) = 0$
-   [Therefore,
    $E(Y(s_i) Y(s_j)) = E(Y(s_i) Y(s_j)) - E(Y(s_i)) E(Y(s_j)) = Cov(Y(s_i), Y(s_j))$]{.orange}
-   We have not defined the function $C(\cdot, \cdot)$ yet
-   $C(\cdot, \cdot)$ tells us *all we need to know* about how our
    spatial random variables are dependent with each other.
    [Why?]{.orange}

## Building Gaussian Processes

-   $C(\cdot, \cdot)$ tells us *all we need to know* about how our
    spatial random variables are dependent with each other.
    [Why?]{.orange}
-   Because we have made a Gaussian assumption
-   Covariance matrix completely defines a mean-zero MVN
-   Covariance function tells us how to build a covariance matrix on any
    set of locations
-   **GP modeling boils down to covariance modeling**
-   Our assumptions on $C(\cdot, \cdot)$ will direct our inference
-   Typically, $C(\cdot, \cdot)$ is a *parametric* function of a small
    number of unknown parameters
-   Covariance parameters describe the underlying process's variance,
       smoothness, spatial range
       -   Let $\theta$ be a vector of covariance parameters, then we can write
       our covariance function as $C_{\theta}(\cdot, \cdot)$
         
         
         ## Covariance function vs Covariance between random variables
         
         - Note we are constructing a covariance function that takes *pairs of spatial locations* as inputs
       - We use that function to model the *covariance between random variables at that pair of spatial locations* 
         - In other words we are making this assumption:
         
         $$Cov(Y(s), Y(s')) = f(s, s')$$
                 for some function $f$ of our choice. In other words we model covariance *uniquely* via the random variables' spatial locations.



## Exponential covariance model

-   Let's consider a simple covariance model
               
               $$ C(s, s') = \sigma^2 \exp \{ -\phi \| s-s' \| \} \qquad \text{Exponential covariance} $$ 
  
  - Assume $Y(\cdot) \sim GP(0, C(\cdot, \cdot))$ with $C(\cdot, \cdot)$ as
defined above 
- On the set of locations
$\cal T = \{\ell_1, \dots, \ell_n\}$ we have

$$ Y_{\cal T} \sim MVN(0, C_{\cal T})  $$ where the $i,j$th element of
$C_{\cal T}$ is
$C_{\cal T}(i,j) = C(\ell_i, \ell_j) = \sigma^2 \exp\{ -\phi \| s_i - s_j \| \}$
  
  
  
  ## Exponential covariance model
  
  ```{r}
set.seed(696)
nobs <- 2000
coords <- data.frame(xcoord = runif(nobs, 0, 1),
                     ycoord = runif(nobs, 0, 1))

head(coords, 5)
```

For example, distance between first and third location,
$$\| s_1 - s_3 \| = \sqrt{(s_{1,x}-s_{3,x})^2 + (s_{1,y}-s_{3,y})^2}$$
  
  ```{r}
sqrt(sum( (coords[1,]-coords[3,])^2 )) 
```

Calculate all pairwise distances

```{r}
Dmat <- as.matrix(dist(coords))
Dmat[1,3]
```

## Exponential covariance model

Fix covariance parameters for simulation

```{r}
sigma2 <- 4
phi <- 3
```

Covariance between $Y(s_1)$ and $Y(s_3)$
  
  ```{r}
sigma2 * exp(-phi * Dmat[1,3])
```

Compute covariance element-wise for all locations

```{r}
C <- sigma2 * exp(- phi * Dmat)
C[1:5, 1:5]
```

## Exponential covariance model

Compute Cholesky decomposition of covariance matrix `C`

```{r}
L <- t(chol(C)) # chol() computes upper cholesky U=t(L), C = crossprod(U) = tcrossprod(L)
```

Sample the MVN

```{r}
Y <- L %*% rnorm(nobs) + 0.1*rnorm(nobs)
```

Plot the data

```{r echo=FALSE}
df <- data.frame(coords, Y)

library(tidyverse)
library(scico)
library(gridExtra)

data_plot <- ggplot(df, aes(xcoord, ycoord, color=Y)) + 
  geom_point(size=1) +
  scale_color_scico(palette="bamO") +
  theme_minimal()

sv <- df %>% with(geoR::variog(data = Y, coords = cbind(xcoord, ycoord), messages=FALSE))
sv_df <- data.frame(dists = sv$u, variogram = sv$v, npairs = sv$n, sd = sv$sd)
sv_plot <- ggplot(sv_df, aes(x=dists, y=variogram)) + geom_point(size=1, shape=8) +
  theme_minimal()

grid.arrange(data_plot, sv_plot, nrow=1)
```


## Exponential covariance model

Now we can keep sampling from *the same Gaussian Process*
  
  ```{r}
Y <- L %*% rnorm(nobs)
```

Plot the data

```{r echo=FALSE}
df <- data.frame(coords, Y)

data_plot <- ggplot(df, aes(xcoord, ycoord, color=Y)) + 
  geom_point(size=1) +
  scale_color_scico(palette="bamO") +
  theme_minimal()

sv <- df %>% with(geoR::variog(data = Y, coords = cbind(xcoord, ycoord), messages=FALSE))
sv_df <- data.frame(dists = sv$u, variogram = sv$v, npairs = sv$n, sd = sv$sd)
sv_plot <- ggplot(sv_df, aes(x=dists, y=variogram)) + geom_point(size=1, shape=8) +
  theme_minimal()

grid.arrange(data_plot, sv_plot, nrow=1)
```

## Exponential covariance model

Now we can keep sampling from *the same Gaussian Process* (at the same locations)

```{r}
Y <- L %*% rnorm(nobs)
```

Plot the data

```{r echo=FALSE}
df <- data.frame(coords, Y)

data_plot <- ggplot(df, aes(xcoord, ycoord, color=Y)) + 
  geom_point(size=1) +
  scale_color_scico(palette="bamO") +
  theme_minimal()

sv <- df %>% with(geoR::variog(data = Y, coords = cbind(xcoord, ycoord), messages=FALSE))
sv_df <- data.frame(dists = sv$u, variogram = sv$v, npairs = sv$n, sd = sv$sd)
sv_plot <- ggplot(sv_df, aes(x=dists, y=variogram)) + geom_point(size=1, shape=8) +
  theme_minimal()

grid.arrange(data_plot, sv_plot, nrow=1)
```

## Exponential covariance model

Now we can keep sampling from *the same Gaussian Process* (at the same locations)

```{r}
Y <- L %*% rnorm(nobs)
```

Plot the data

```{r echo=FALSE}
df <- data.frame(coords, Y)

data_plot <- ggplot(df, aes(xcoord, ycoord, color=Y)) + 
  geom_point(size=1) +
  scale_color_scico(palette="bamO") +
  theme_minimal()

sv <- df %>% with(geoR::variog(data = Y, coords = cbind(xcoord, ycoord), messages=FALSE))
sv_df <- data.frame(dists = sv$u, variogram = sv$v, npairs = sv$n, sd = sv$sd)
sv_plot <- ggplot(sv_df, aes(x=dists, y=variogram)) + geom_point(size=1, shape=8) +
  theme_minimal()

grid.arrange(data_plot, sv_plot, nrow=1)
```

## Exponential covariance model

Now we can keep sampling from *the same Gaussian Process* (at the same locations)

```{r}
Y <- L %*% rnorm(nobs)
```

Plot the data

```{r echo=FALSE}
df <- data.frame(coords, Y)

data_plot <- ggplot(df, aes(xcoord, ycoord, color=Y)) + 
  geom_point(size=1) +
  scale_color_scico(palette="bamO") +
  theme_minimal()

sv <- df %>% with(geoR::variog(data = Y, coords = cbind(xcoord, ycoord), messages=FALSE))
sv_df <- data.frame(dists = sv$u, variogram = sv$v, npairs = sv$n, sd = sv$sd)
sv_plot <- ggplot(sv_df, aes(x=dists, y=variogram)) + geom_point(size=1, shape=8) +
  theme_minimal()

grid.arrange(data_plot, sv_plot, nrow=1)
```

## Exponential covariance model

Now we can keep sampling from *the same Gaussian Process* (at the same locations)

```{r}
Y <- L %*% rnorm(nobs)
```

Plot the data

```{r echo=FALSE}
df <- data.frame(coords, Y)

data_plot <- ggplot(df, aes(xcoord, ycoord, color=Y)) + 
  geom_point(size=1) +
  scale_color_scico(palette="bamO") +
  theme_minimal()

sv <- df %>% with(geoR::variog(data = Y, coords = cbind(xcoord, ycoord), messages=FALSE))
sv_df <- data.frame(dists = sv$u, variogram = sv$v, npairs = sv$n, sd = sv$sd)
sv_plot <- ggplot(sv_df, aes(x=dists, y=variogram)) + geom_point(size=1, shape=8) +
  theme_minimal()

grid.arrange(data_plot, sv_plot, nrow=1)
```

## Exponential covariance model

Now we can keep sampling from *the same Gaussian Process* (at the same locations)

```{r}
Y <- L %*% rnorm(nobs)
```

Plot the data

```{r echo=FALSE}
df <- data.frame(coords, Y)

data_plot <- ggplot(df, aes(xcoord, ycoord, color=Y)) + 
  geom_point(size=1) +
  scale_color_scico(palette="bamO") +
  theme_minimal()

sv <- df %>% with(geoR::variog(data = Y, coords = cbind(xcoord, ycoord), messages=FALSE))
sv_df <- data.frame(dists = sv$u, variogram = sv$v, npairs = sv$n, sd = sv$sd)
sv_plot <- ggplot(sv_df, aes(x=dists, y=variogram)) + geom_point(size=1, shape=8) +
  theme_minimal()

grid.arrange(data_plot, sv_plot, nrow=1)
```


## Stationary covariance functions

-  $C(\cdot, \cdot)$ is **stationary** if $C(s, s + h) = C(h)$ for all $h \in \Re^d$ such that $s, s+h \in D$
  -  In other words, covariance only depends on the separation vector $h$
  -  In other words, $C(s, s') = C(s-s')$, we are modeling covariance as not depending on the *absolute* location of $s$ and $s'$ but only their *relative* position in the spatial domain


## Isotropic covariance functions

-  $C(\cdot, \cdot)$ is **isotropic** if $C(s, s + h) = C(||h||)$ for all $h \in \Re^d$ such that $s, s+h \in D$
-  In other words, covariance only depends on the *length* of the separation vector $h$
-  In other words, $C(s, s') = C(||s-s'||)$, we are modeling covariance as not depending on the *absolute* location of $s$ and $s'$ but only their *relative* position in the spatial domain, without regard to orientation
                                 -  We can write $C(s, s')=C(d)$ where $d = ||s-s'||$ is the distance between $s$ and $s'$


## Example: Exponential covariance

```{r echo=FALSE}
library(scico)
library(gridExtra)
library(magrittr)

coords <- expand.grid(gg <- seq(0,1,length.out=50), gg)
m_coords <- as.matrix(coords)
ch <- coords %>% filter(Var1 == 0)
```

$$ C(d) = \sigma^2 \exp \{ -\phi d \} $$ 

 - $\phi$ is the *spatial decay* parameter. Sometimes we use $\alpha=1/\phi$, called the *range* 
 - $\sigma^2$ is the spatial variance or *sill*
 
Example with $\sigma^2=2,\ \phi=10$:

```{r echo=FALSE}
sigmasq <- 2; phi <- 10; nu <- 0.5; tausq <- 0

exponential_covar <- gg %>% sapply(\(x) meshed:::Cov_matern_h(x, sigmasq, phi, nu, tausq))

cov_df <- data.frame(h = gg, Exponential=exponential_covar)
cov_plot <- ggplot(cov_df, aes(h, Exponential)) +
  geom_line() +
  theme_minimal() +
  labs(x="Distance", y="Covariance")

# manually do everything
#maxmin_order <- GpGp::order_maxmin(m_coords)
#mo_coords <- m_coords[maxmin_order,]
#nneighbors <- GpGp::find_ordered_nn(mo_coords, m=20)
#test <- GpGp::vecchia_Linv(c(sigmasq, 1/phi, nu, tausq), "matern_isotropic", mo_coords, nneighbors)

set.seed(1)
w <- GpGp::fast_Gp_sim(c(sigmasq, 1/phi, nu, tausq), "matern_isotropic", m_coords, m = 20)
df <- data.frame(coords, w)

data_plot <- ggplot(df, aes(Var1, Var2, fill=w)) + 
  geom_raster() +
  scale_fill_scico(palette="lapaz") +
  labs(fill="Y", x="x-coord", y="y-coord") + theme_minimal()

grid.arrange(data_plot, cov_plot, nrow=1)
```

## Example: Squared exponential covariance (aka Gaussian covariance)

$$ C(d) = \sigma^2 \exp \{ -\phi d^2 \} $$


- $\phi$ is the *spatial decay* parameter. Sometimes we use $\alpha=1/\phi$, called the *range* 
- $\sigma^2$ is the spatial variance or *sill*
 
Example with $\sigma^2=2,\ \phi=20$:
 
```{r echo=FALSE}
sigmasq <- 2; phi <- 20; nu <- 4.5; tausq <- 0

sqexp_covar <- gg %>% sapply(\(x) sigmasq * exp(-phi * x^2))

cov_df %<>% mutate(Sq.Exponential=sqexp_covar)

cov_plot <- ggplot(cov_df, aes(h, Sq.Exponential)) +
  geom_line() +
  theme_minimal() +
  labs(x="Distance", y="Covariance")

set.seed(1)
# cheating a little but it would look pretty much the same ;)
w <- GpGp::fast_Gp_sim(c(sigmasq, 1/phi, 4.5, tausq), "matern_isotropic", m_coords, m = 20)
df <- data.frame(coords, w)

data_plot <- ggplot(df, aes(Var1, Var2, fill=w)) + 
  geom_raster() +
  scale_fill_scico(palette="lapaz") +
  labs(fill="Y", x="x-coord", y="y-coord") + theme_minimal()

grid.arrange(data_plot, cov_plot, nrow=1)
```


## Example: the Matern model of covariance

$$ C(d) = \sigma^2 \frac{2^{1-\nu}}{\Gamma(\nu)} (\phi d)^{\nu}   K_{\nu}( \phi d^2 )$$

- $\Gamma(\cdot)$ is the Gamma function, $K_{\nu}(\cdot)$ is a modified Bessel function
- $\nu$ is the smoothness parameter
- If $\nu=0.5$ this is the exponential covariance model
- At $\nu\to\infty$ this becomes the squared exponential covariance model

## Example: the Matern model of covariance

$$ C(d) = \sigma^2 \frac{2^{1-\nu}}{\Gamma(\nu)} (\phi d)^{\nu}   K_{\nu}( \phi d^2 )$$
Example with $\sigma^2 = 2,\ \phi=15,\ \nu=1.2$:

```{r echo=FALSE}
sigmasq <- 2; phi <- 15; nu <- 1.2; tausq <- 0

matern_covar <- gg %>% sapply(\(x) meshed:::Cov_matern_h(x, sigmasq, phi, nu, tausq))

cov_df %<>% mutate(Matern12=matern_covar)

cov_plot <- ggplot(cov_df, aes(h, Matern12)) +
  geom_line() +
  theme_minimal() +
  labs(x="Distance", y="Covariance")

set.seed(1)
# cheating a little but it would look pretty much the same ;)
w <- GpGp::fast_Gp_sim(c(sigmasq, 1/phi, 4.5, tausq), "matern_isotropic", m_coords, m = 20)
df <- data.frame(coords, w)

data_plot <- ggplot(df, aes(Var1, Var2, fill=w)) + 
  geom_raster() +
  scale_fill_scico(palette="lapaz") +
  labs(fill="Y", x="x-coord", y="y-coord") + theme_minimal()

grid.arrange(data_plot, cov_plot, nrow=1)
```

```{r echo=FALSE}
library(rootSolve)
sigmasq <- 1; tausq <- 0

nu_vals <- c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, Inf)
nu_colors <- colorRampPalette(c("royalblue", "orange"))(length(nu_vals))

phi_target <- function(nu){
  vecmatern <- function(x) sapply(x, \(s) ifelse(nu==Inf, 
                                                 sigmasq * exp(-s * 0.25^2), 
                                                 meshed:::Cov_matern_h(0.25, sigmasq, s, nu, tausq)))
  mymatern <- function(x) vecmatern(x)-0.5
  target_phi <- rootSolve::uniroot.all(mymatern, lower=1, upper=50)
return(target_phi)
}

make_plot <- function(nuseq){
  
  matern_covar <- bind_rows(nuseq %>% lapply(\(nu) {
  cova <- gg %>% sapply(\(x) ifelse(nu == Inf,
                                    sigmasq * exp(-phi_target(nu) * x^2),
                                    meshed:::Cov_matern_h(x, sigmasq, phi_target(nu), nu, tausq)))
  return(data.frame(Distance=gg, nu=nu, Covariance=cova))
  }))

  cov_plot <- ggplot(matern_covar, aes(x=Distance, y=Covariance, color=factor(nu))) +
    geom_line() +
    theme_minimal() +
    scale_discrete_manual(aesthetics = "color", values=nu_colors[1:length(nuseq)]) +
    coord_cartesian(ylim=c(0,1)) +
    labs(x="Distance", y="Covariance", color="Smoothness")
  
  return(cov_plot)
}
```


## Comparing Matern covariances with different smoothness

```{r echo=FALSE}
make_plot(nu_vals[1])
```


## Comparing Matern covariances with different smoothness

```{r echo=FALSE}
make_plot(nu_vals[1:3])
```


## Comparing Matern covariances with different smoothness

```{r echo=FALSE}
make_plot(nu_vals[1:6])
```


## Comparing Matern covariances with different smoothness

```{r echo=FALSE}
make_plot(nu_vals[1:11])
```


## Comparing Matern covariances with different smoothness

```{r echo=FALSE}
make_plot(nu_vals) 
```




## Other examples

-  Power exponential covariance
$$ C(d) = \sigma^2 \exp \{ -\phi d^\gamma \}  \qquad \text{where}\qquad \gamma\in (0,2] $$ 
-  Exponential covariance with nugget term
$$ C(d) = \sigma^2 \exp\{ -\phi d \} + \cal I_{d=0}(\tau^2)$$
    - $\tau^2$ is called the *nugget*, whereas $\sigma^2$ is the *sill*
    - the nugget term represents measurement error
    - we can add a nugget term to any covariance function

-  See textbook at page 28 for other examples


## Nugget covariance

```{r echo=FALSE}
sigmasq <- 1; phi <- 5; nu <- 1.2; tausq <- 0

matern_covar <- gg %>% sapply(\(x) meshed:::Cov_matern_h(x, sigmasq, phi, nu, tausq))

cov_df <- data.frame(h=gg, Covariance=matern_covar)

(cov_plot <- ggplot(cov_df, aes(h, Covariance)) +
  geom_line() +
  theme_minimal() +
  labs(x="Distance", y="Covariance") +
  coord_cartesian(ylim=c(0, 1.5)) +
  geom_point(data=data.frame(h=0, Covariance=1.2), color="black"))
```

## Stationary *anisotropic* covariance functions

- All covariance functions need to be symmetric, so $C(s, s') = C(s', s)$
- Stationary *isotropic* covariance function $C(s, s') = C(\|h\|)$ depends only on distance $\|h \| = \sqrt{ (s - s')^T (s- s') }$
                                   - Stationary *anisotropic*? Just need $C(s, s') = C(s-s') = C(h)$
                                   - Suppose $A$ is a $d\times d$ dimensional positive definite matrix 
                                 - Simple way to introduce anisotropy: make the covariance function depend on 
                                 
                                 $$ \sqrt{ (s - s')^T A (s- s') } $$
                                   
                                   - Obtain the isotropic case by letting $A = \phi I_d$, with $\phi > 0$.
                                 - Example: anisotropic exponential covariance
                                 
                                 $$ C(h) = \sigma^2 \exp \left\{ - \sqrt{h^T A h} \right\} $$
                                   
                                   - Example: anisotropic power exponential covariance
                                 
                                 $$ C(h) = \sigma^2 \exp \left\{ - (h^T A h)^{\frac{\gamma}{2}}  \right\} $$
                                   
                                   ## Example: anisotropic power exponential covariance
                                   
                                   -  Let's simulate data with this covariance function:

$$ C(h) = \sigma^2 \exp \left\{ - (h^T A h)^{\frac{\gamma}{2}} \right\} $$

-  We have $A = \begin{bmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{bmatrix}$
-  If we let $a_{ij} \in \Re$ without constraint, we might break symmetry and positive definiteness of $A$ which we require
-  [How could we define our covariance function so that the resulting $A$ matrix is always symmetric positive definite?]{.orange}
-  [What we need to do is *reparametrize* our covariance function]{.orange}

## Example: anisotropic power exponential covariance

-  [How could we define our covariance function so that the resulting $A$ matrix is always symmetric positive definite?]{.orange}
-  [What we need to do is *reparametrize* our covariance function]{.orange}
-  One option: Cholesky. $L_A = \begin{bmatrix} l_{11} & 0 \\ l_{21} & l_{22} \end{bmatrix}$, then 

$$ C(h) = \sigma^2 \exp \left\{ - (h^T L_A L_A^T h )^{\frac{\gamma}{2}} \right\} $$

-  The matrix $L_A L_A^T$ is symmetric positive definite by construction
-  If $l_{21} = l_{22} = 1$ then $L_A L_A^T = I_d$ and we are back to the isotropic case
-  But this is not very interpretable. How do we assign a "meaning" to $l_{11}$ and $l_{22}$? 
-  Instead, define $S = \text{diag}\{ \sqrt{\phi}, \sqrt{s \phi} \}$, $R = \begin{bmatrix} 1 & \rho \\ \rho & 1 \end{bmatrix}$, where $s, \phi > 0$, $-1<\rho<1$. Then
$$ C(h) = \sigma^2 \exp \left\{ - (h^T S R S h )^{\frac{\gamma}{2}} \right\} $$



```{r echo=FALSE}
coords <- expand.grid(xgrid <- seq(0, 1, length.out=114), ygrid <- seq(0, 1, length.out=57))
colnames(coords) <- c("xcoord", "ycoord")
cmat <- coords %>% as.matrix()

mycppcode <- "
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace std;

//[[Rcpp::export]]
arma::mat anisotropic_powerexp(const arma::mat& cx, const arma::mat& cy, const arma::vec& param){

  double sigmasq = param(0);
  double phi = param(1);
  double s = param(2);
  double rho = param(3);
  
  arma::mat S = arma::eye(2,2);
  S(0,0) = sqrt(phi);
  S(1,1) = sqrt(s);
  arma::mat R = arma::eye(2,2);
  R(1,0) = rho;
  R(0,1) = rho;
  arma::mat A = S * R * S;
  
  double expo = param(4);
  double nugget = param(5);
  
  arma::mat result = arma::zeros(cx.n_rows, cy.n_rows);
  for(unsigned int i=0; i<cx.n_rows; i++){
    for(unsigned int j=0; j<cy.n_rows; j++){
    
      arma::rowvec h = cx.row(i) - cy.row(j);
      double inner = arma::conv_to<double>::from( sqrt(h * A * h.t()) );
      result(i, j) = sigmasq * exp(- pow( inner, expo ) );
      if(arma::all(h == 0)){
        result(i, j) += nugget;
      }
    }
  }
  return result;
}
"
Rcpp::sourceCpp(code=mycppcode)

make_anisotropic_plot <- function(params){
  C <- anisotropic_powerexp(cmat, cmat, params)
  L <- t(chol(C))
  set.seed(1)
  Y <- L %*% rnorm(nrow(coords))

  df <- data.frame(coords, Y)
  data_plot <- ggplot(df, aes(xcoord, ycoord, fill=Y)) + 
    geom_raster() +
    scale_fill_scico(palette="bamO") +
    theme_void()
  return(data_plot)
}
```

## Example: anisotropic power exponential covariance

```{r echo=FALSE}
# sigmasq, phi, s, rho, exponent, nugget
params <- c(1, 1, 1, 0, 1, 1e-5)
```

$C(h) = \sigma^2 \exp \left\{ - \left(h^T \begin{bmatrix} \sqrt{\phi} & 0 \\ 0 &  \sqrt{s \phi} \end{bmatrix} \begin{bmatrix} 1 & \rho \\ \rho & 1 \end{bmatrix} \begin{bmatrix} \sqrt{\phi} & 0 \\ 0 &  \sqrt{s \phi} \end{bmatrix} h \right)^{\frac{\gamma}{2}} \right\}$, with $$\sigma^2 = `r params[1]`,\ \phi = `r params[2]`,\ s = `r params[3]`,\ \rho = `r params[4]`,\ \gamma = `r params[5]`$$

```{r echo=FALSE, fig.asp=0.5}
make_anisotropic_plot(params)
```


## Example: anisotropic power exponential covariance

```{r echo=FALSE}
# sigmasq, phi, s, rho, exponent, nugget
params <- c(1, 5, 1, 0, 1.5, 1e-5)
```

$C(h) = \sigma^2 \exp \left\{ - \left(h^T \begin{bmatrix} \sqrt{\phi} & 0 \\ 0 &  \sqrt{s \phi} \end{bmatrix} \begin{bmatrix} 1 & \rho \\ \rho & 1 \end{bmatrix} \begin{bmatrix} \sqrt{\phi} & 0 \\ 0 &  \sqrt{s \phi} \end{bmatrix} h \right)^{\frac{\gamma}{2}} \right\}$, with $$\sigma^2 = `r params[1]`,\ \phi = `r params[2]`,\ s = `r params[3]`,\ \rho = `r params[4]`,\ \gamma = `r params[5]`$$

```{r echo=FALSE, fig.asp=0.5}
make_anisotropic_plot(params)
```

## Example: anisotropic power exponential covariance

```{r echo=FALSE}
# sigmasq, phi, s, rho, exponent, nugget
params <- c(1, 5, 3, 0, 1.5, 1e-5)
```

$C(h) = \sigma^2 \exp \left\{ - \left(h^T \begin{bmatrix} \sqrt{\phi} & 0 \\ 0 &  \sqrt{s \phi} \end{bmatrix} \begin{bmatrix} 1 & \rho \\ \rho & 1 \end{bmatrix} \begin{bmatrix} \sqrt{\phi} & 0 \\ 0 &  \sqrt{s \phi} \end{bmatrix} h \right)^{\frac{\gamma}{2}} \right\}$, with $$\sigma^2 = `r params[1]`,\ \phi = `r params[2]`,\ s = `r params[3]`,\ \rho = `r params[4]`,\ \gamma = `r params[5]`$$

```{r echo=FALSE, fig.asp=0.5}
make_anisotropic_plot(params)
```


## Example: anisotropic power exponential covariance

```{r echo=FALSE}
# sigmasq, phi, s, rho, exponent, nugget
params <- c(1, 5, 10, 0, 1.5, 1e-5)
```

$C(h) = \sigma^2 \exp \left\{ - \left(h^T \begin{bmatrix} \sqrt{\phi} & 0 \\ 0 &  \sqrt{s \phi} \end{bmatrix} \begin{bmatrix} 1 & \rho \\ \rho & 1 \end{bmatrix} \begin{bmatrix} \sqrt{\phi} & 0 \\ 0 &  \sqrt{s \phi} \end{bmatrix} h \right)^{\frac{\gamma}{2}} \right\}$, with $$\sigma^2 = `r params[1]`,\ \phi = `r params[2]`,\ s = `r params[3]`,\ \rho = `r params[4]`,\ \gamma = `r params[5]`$$

```{r echo=FALSE, fig.asp=0.5}
make_anisotropic_plot(params)
```



## Example: anisotropic power exponential covariance

```{r echo=FALSE}
# sigmasq, phi, s, rho, exponent, nugget
params <- c(1, 5, 25, 0, 1.5, 1e-5)
```

$C(h) = \sigma^2 \exp \left\{ - \left(h^T \begin{bmatrix} \sqrt{\phi} & 0 \\ 0 &  \sqrt{s \phi} \end{bmatrix} \begin{bmatrix} 1 & \rho \\ \rho & 1 \end{bmatrix} \begin{bmatrix} \sqrt{\phi} & 0 \\ 0 &  \sqrt{s \phi} \end{bmatrix} h \right)^{\frac{\gamma}{2}} \right\}$, with $$\sigma^2 = `r params[1]`,\ \phi = `r params[2]`,\ s = `r params[3]`,\ \rho = `r params[4]`,\ \gamma = `r params[5]`$$

```{r echo=FALSE, fig.asp=0.5}
make_anisotropic_plot(params)
```

## Example: anisotropic power exponential covariance

```{r echo=FALSE, fig.asp=0.5}
# sigmasq, phi, s, rho, exponent, nugget
params <- c(1, 5, 25, 0.3, 1.5, 1e-5)
```

$C(h) = \sigma^2 \exp \left\{ - \left(h^T \begin{bmatrix} \sqrt{\phi} & 0 \\ 0 &  \sqrt{s \phi} \end{bmatrix} \begin{bmatrix} 1 & \rho \\ \rho & 1 \end{bmatrix} \begin{bmatrix} \sqrt{\phi} & 0 \\ 0 &  \sqrt{s \phi} \end{bmatrix} h \right)^{\frac{\gamma}{2}} \right\}$, with $$\sigma^2 = `r params[1]`,\ \phi = `r params[2]`,\ s = `r params[3]`,\ \rho = `r params[4]`,\ \gamma = `r params[5]`$$

```{r echo=FALSE, fig.asp=0.5}
make_anisotropic_plot(params)
```

## Example: anisotropic power exponential covariance

```{r echo=FALSE}
# sigmasq, phi, s, rho, exponent, nugget
params <- c(1, 5, 25, 0.9, 1.5, 1e-5)
```

$C(h) = \sigma^2 \exp \left\{ - \left(h^T \begin{bmatrix} \sqrt{\phi} & 0 \\ 0 &  \sqrt{s \phi} \end{bmatrix} \begin{bmatrix} 1 & \rho \\ \rho & 1 \end{bmatrix} \begin{bmatrix} \sqrt{\phi} & 0 \\ 0 &  \sqrt{s \phi} \end{bmatrix} h \right)^{\frac{\gamma}{2}} \right\}$, with $$\sigma^2 = `r params[1]`,\ \phi = `r params[2]`,\ s = `r params[3]`,\ \rho = `r params[4]`,\ \gamma = `r params[5]`$$

```{r echo=FALSE, fig.asp=0.5}
make_anisotropic_plot(params)
```


## Example: anisotropic power exponential covariance

```{r echo=FALSE}
params <- c(1, 5, 25, -.9, 1.5, 1e-5)
```

$C(h) = \sigma^2 \exp \left\{ - \left(h^T \begin{bmatrix} \sqrt{\phi} & 0 \\ 0 &  \sqrt{s \phi} \end{bmatrix} \begin{bmatrix} 1 & \rho \\ \rho & 1 \end{bmatrix} \begin{bmatrix} \sqrt{\phi} & 0 \\ 0 &  \sqrt{s \phi} \end{bmatrix} h \right)^{\frac{\gamma}{2}} \right\}$, with $$\sigma^2 = `r params[1]`,\ \phi = `r params[2]`,\ s = `r params[3]`,\ \rho = `r params[4]`,\ \gamma = `r params[5]`$$

```{r echo=FALSE, fig.asp=0.5}
make_anisotropic_plot(params)
```


## Visualizing a stationary anisotropic covariance function

```{r echo=FALSE}
(cov_plot <- ggplot(cov_df, aes(h, Covariance)) +
  geom_line() +
  theme_minimal() +
  labs(x="Distance", y="Covariance"))
```

- In the stationary anisotropic cases we considered earlier, $C(s, s')$ depended on the magnitude of $s-s'$ (ie the distance between $s$ and $s'$) as well as the orientation
       - Therefore there is at least one new dimension to the plot
       
       ```{r echo=FALSE}
       mycppcode <- "
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace std;

//[[Rcpp::export]]
double anisotropic_powerexp_polarcoords(double r, double degrees, const arma::vec& param){

  double sigmasq = param(0);
  double phi = param(1);
  double s = param(2);
  double rho = param(3);
  
  double angle = degrees * M_PI/180.0;
  
  arma::mat S = arma::eye(2,2);
  S(0,0) = sqrt(phi);
  S(1,1) = sqrt(s);
  arma::mat R = arma::eye(2,2);
  R(1,0) = rho;
  R(0,1) = rho;
  arma::mat A = S * R * S;
  
  double expo = param(4);
  double nugget = param(5);
  
  arma::vec h = arma::zeros(2);
  h(0) = r * cos(angle);
  h(1) = r * sin(angle);
  
  double inner = arma::conv_to<double>::from( sqrt(h.t() * A * h) );
  double covariance = sigmasq * exp(- pow( inner, expo ) );
  if(arma::all(h == 0)){
    covariance += nugget;
  }

  return covariance;
}
"
Rcpp::sourceCpp(code=mycppcode)
```

## Visualizing a stationary anisotropic covariance function

```{r echo=FALSE}
# sigmasq, phi, s, rho, exponent, nugget
params <- c(1, 5, 25, 0.3, 1.5, 1e-5)
```

$C(h) = \sigma^2 \exp \left\{ - \left(h^T \begin{bmatrix} \sqrt{\phi} & 0 \\ 0 &  \sqrt{s \phi} \end{bmatrix} \begin{bmatrix} 1 & \rho \\ \rho & 1 \end{bmatrix} \begin{bmatrix} \sqrt{\phi} & 0 \\ 0 &  \sqrt{s \phi} \end{bmatrix} h \right)^{\frac{\gamma}{2}} \right\}$, with $$\sigma^2 = `r params[1]`,\ \phi = `r params[2]`,\ s = `r params[3]`,\ \rho = `r params[4]`,\ \gamma = `r params[5]`$$
  
  ```{r echo=FALSE}
coords <- expand.grid(xgrid <- seq(0, 1, length.out=80), ygrid <- seq(0, 1, length.out=80))
colnames(coords) <- c("xcoord", "ycoord")
cmat <- coords %>% as.matrix()

data_plot <- make_anisotropic_plot(params)


polar_grids <- list(d = seq(0, sqrt(2), length.out=200), 
                    degrees = seq(0, 360, 14))
polar_df <- expand.grid(polar_grids)

ani_cov_df <- polar_df %>% 
  rowwise() %>% 
  mutate(covariance = anisotropic_powerexp_polarcoords(d, degrees, params))

degree_colors <- colorRampPalette(c("royalblue", "orange"))(length(polar_grids$degrees))


ani_cov_plot <- ggplot(ani_cov_df, aes(d, covariance, color=factor(degrees))) +
  geom_line() +
  theme_minimal() +
  scale_discrete_manual(aesthetics = "color", values=degree_colors) +
  labs(x="Distance", y="Covariance", color="Degrees") +
  coord_cartesian(ylim=c(0, 1.1))

grid.arrange(data_plot, ani_cov_plot, nrow=1)
```


## Visualizing a stationary anisotropic covariance function

```{r echo=FALSE}
# sigmasq, phi, s, rho, exponent, nugget
params <- c(1, 5, 25, 0.3, 1.5, 1e-5)
```

$C(h) = \sigma^2 \exp \left\{ - \left(h^T \begin{bmatrix} \sqrt{\phi} & 0 \\ 0 &  \sqrt{s \phi} \end{bmatrix} \begin{bmatrix} 1 & \rho \\ \rho & 1 \end{bmatrix} \begin{bmatrix} \sqrt{\phi} & 0 \\ 0 &  \sqrt{s \phi} \end{bmatrix} h \right)^{\frac{\gamma}{2}} \right\}$, with $$\sigma^2 = `r params[1]`,\ \phi = `r params[2]`,\ s = `r params[3]`,\ \rho = `r params[4]`,\ \gamma = `r params[5]`$$
  
  ```{r echo=FALSE}
coords <- expand.grid(xgrid <- seq(0, 1, length.out=80), ygrid <- seq(0, 1, length.out=80))
colnames(coords) <- c("xcoord", "ycoord")
cmat <- coords %>% as.matrix()

data_plot <- make_anisotropic_plot(params)


polar_grids <- list(d = seq(0, sqrt(2), length.out=10), 
                    degrees = seq(0, 360, length.out=200))
polar_df <- expand.grid(polar_grids)

ani_cov_df <- polar_df %>% 
  rowwise() %>% 
  mutate(covariance = anisotropic_powerexp_polarcoords(d, degrees, params))

d_colors <- colorRampPalette(c("royalblue", "orange"))(length(polar_grids$d))


ani_cov_plot <- ggplot(ani_cov_df, aes(degrees, covariance, color=factor(round(d,2)))) +
  geom_line() +
  theme_minimal() +
  scale_discrete_manual(aesthetics = "color", values=d_colors) +
  labs(x="Degrees", y="Covariance", color="Distance") +
  coord_cartesian(ylim=c(0, 1.1))

grid.arrange(data_plot, ani_cov_plot, nrow=1)
```

## Visualizing a stationary anisotropic covariance function

```{r echo=FALSE}
# sigmasq, phi, s, rho, exponent, nugget
params <- c(1, 5, 5, 0.9, 1.5, 1e-5)
```

$C(h) = \sigma^2 \exp \left\{ - \left(h^T \begin{bmatrix} \sqrt{\phi} & 0 \\ 0 &  \sqrt{s \phi} \end{bmatrix} \begin{bmatrix} 1 & \rho \\ \rho & 1 \end{bmatrix} \begin{bmatrix} \sqrt{\phi} & 0 \\ 0 &  \sqrt{s \phi} \end{bmatrix} h \right)^{\frac{\gamma}{2}} \right\}$, with $$\sigma^2 = `r params[1]`,\ \phi = `r params[2]`,\ s = `r params[3]`,\ \rho = `r params[4]`,\ \gamma = `r params[5]`$$
  
  ```{r echo=FALSE}
coords <- expand.grid(xgrid <- seq(0, 1, length.out=80), ygrid <- seq(0, 1, length.out=80))
colnames(coords) <- c("xcoord", "ycoord")
cmat <- coords %>% as.matrix()

data_plot <- make_anisotropic_plot(params)


polar_grids <- list(d = seq(0, sqrt(2), length.out=10), 
                    degrees = seq(0, 360, length.out=200))
polar_df <- expand.grid(polar_grids)

ani_cov_df <- polar_df %>% 
  rowwise() %>% 
  mutate(covariance = anisotropic_powerexp_polarcoords(d, degrees, params))

d_colors <- colorRampPalette(c("royalblue", "orange"))(length(polar_grids$d))


ani_cov_plot <- ggplot(ani_cov_df, aes(degrees, covariance, color=factor(round(d,2)))) +
  geom_line() +
  theme_minimal() +
  scale_discrete_manual(aesthetics = "color", values=d_colors) +
  labs(x="Degrees", y="Covariance", color="Distance") +
  coord_cartesian(ylim=c(0, 1.1))

grid.arrange(data_plot, ani_cov_plot, nrow=1)
```



## Interpretation of covariance models

-  We have made covariance plot for 
- stationary isotropic covariances
- stationary anisotropic covariances
-  [Discuss ideas for making covariance plots of *nonstationary* covariance models]{.orange}
-  [No specific model in mind, just think of what *nonstationary* implies]{.orange}

## Interpretation of covariance models

-  One reason we choose a specific covariance model is because we think it is a **good model** of spatial variation in the data
-  Another reason we choose a specific covariance model is to be able to **interpret** the results of our estimation
-  There is a **tradeoff** between ease of estimation, interpretation, and flexibility to realistic assumptions
-  Parameters themselves are not always directly interpretable
-  Different covariance functions' parameters may not be directly comparable in how they are interpreted
-  Covariance plot is an effective tool for interpretation
-  Covariance plot becomes more difficult to do the more the covariance function is complicated!



## Interpretation of covariance models

-  Other questions that we could answer:
    - what is the signal-to-noise ratio?
    - how far do we need to take two locations in order for their correlation to be less than 0.95? 
    - what is the distance between locations at which correlation falls under 0.05?
-  The second question is also called the **effective range**

```{r echo=FALSE}
(cov_plot <- ggplot(cov_df, aes(h, Covariance)) +
  geom_line() +
  theme_minimal() +
  labs(x="Distance", y="Covariance") +
  coord_cartesian(ylim=c(0, 1.5)) +
  geom_point(data=data.frame(h=0, Covariance=1.2), color="black"))
```

## Effective range

The distance $d$ such that $Cor(d)=0.05$ where $Cor(\cdot)$ is the implied correlation function.

*Example:* exponential correlation function $Cor(d) = \exp\{ -\phi d \}$.

```{r echo=FALSE}
(cov_plot <- ggplot(cov_df, aes(h, Covariance)) +
  geom_line() +
  theme_minimal() +
  labs(x="Distance", y="Covariance") +
  coord_cartesian(ylim=c(0, 1)) +
  geom_hline(yintercept=0.05, color="red") + 
  geom_vline(xintercept=0.865, lty=2) +
  scale_x_continuous(breaks=c(0, 0.25, 0.50, 0.75, 0.865, 1)) +
  scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.50, 0.75, 1)))
```

[Is this an exponential correlation plot?]{.orange}


## Some ways to obtain other covariance functions

- Sum: if $C_1(s, s')$ and $C_2(s, s')$ are two covariance functions, then 
$$C(s, s') = C_1(s, s') + C_2(s, s')$$ 
  is also a valid covariance function

- Product: if $C_1(s, s')$ and $C_2(s, s')$ are two covariance functions, then 
$$C(s, s') = C_1(s, s') \cdot C_2(s, s')$$ 
is also a valid covariace function


## Some ways to obtain other covariance functions: cont'd
                                  
                                  -  If $w_1(\cdot) \sim GP(0, C_1(\cdot, \cdot))$ is independent of $w_2(\cdot) \sim GP(0, C_2(\cdot, \cdot))$ then $w(\cdot) = w_1(\cdot) + w_2(\cdot) \sim GP(0, C_1(\cdot, \cdot) + C_2(\cdot, \cdot))$
                                    
                                    -  If $w_1(\cdot) \sim GP(0, C_1(\cdot, \cdot))$ is independent of $w_2(\cdot) \sim GP(0, C_2(\cdot, \cdot))$ then $w(\cdot) = w_1(\cdot) \cdot w_2(\cdot)$ is such that $Cov(w(\cdot)) = C(\cdot, \cdot) = C_1(\cdot, \cdot) \cdot C_2(\cdot, \cdot)$. But notice $w(\cdot)$ is not a GP.
                                  
                                  ## Revisiting the nugget term
                                  
                                  $$ C(s, s') = \sigma^2 \exp\{ -\phi \| s-s' \| \} + \cal I_{s=s'}(\tau^2)$$

-  We see that if $w(\cdot) \sim GP(0, C)$ then we can write

$$w(\cdot) = v(\cdot) + \varepsilon(\cdot)$$
where $v(\cdot)$ is a GP with exponential covariance and $\varepsilon(\cdot)$ is a white noise process


## Some ways to obtain other covariance functions: cont'd
                                  
                                  - Direct sum: if $C_1(s_1, s_1')$ is a covariance function on domain $D_1$ and $C_2(s_2, s_2')$ is a covariance function on domain $D_2$ then 
                                  $$C(s, s') = C_1(s_1, s_1') + C_2(s_2, s_2')$$
is a covariance function on the domain $D_1 \times D_2$.

### Example: 

Suppose $D_1=D \subset \Re^2$ is the usual spatial domain, $D_2=T\subset \Re$ is the time domain. 

Let $u(s,t) = w(s) + v(t)$ where $w(\cdot)$ and $v(\cdot)$ are GPs defined on $D$ and $T$, respectively, we have that $u(s,t)$ is a GP on $D \times T$. 

In other words, this is a simple way to define a *spatiotemporal process*

- Tensor product: if $C_1(s_1, s_1')$ is a covariance function on domain $D_1$ and $C_2(s_2, s_2')$ is a covariance function on domain $D_2$ then 
$$C(s, s') = C_1(s_1, s_1') \cdot C_2(s_2, s_2')$$
  is a covariance function on the domain $D_1 \times D_2$
  
  
  ## Examples of non-stationary covariances
  
  -  Paciorek and Schervish (2006), based on the Matern model:
  
  $$C(s_i, s_j) = \sigma^2 \frac{2^{1-\nu}}{\Gamma(\nu)} |\Sigma_i|^{-\frac{1}{4}} |\Sigma_j|^{\frac{1}{4}} \left| \frac{\Sigma_i + \Sigma_j}{2} \right|^{-\frac{1}{2}} \left(2 \sqrt{\nu Q_{ij}} \right)^{\nu} K_{\nu}\left(2 \sqrt{\nu Q_{ij}} \right) $$
    where $Q_{ij} = (s_i - s_j)^T \left(\frac{\Sigma_i + \Sigma_j}{2} \right)^{-1}(s_i - s_j)$ and $\Sigma_i = \Sigma(x_i)$ is called the kernel matrix.
    
    -  Another example 
    
    $$C(s_i, s_j) = \exp\{ \phi(s_i) + \phi(s_j) \} M(s_i, s_j)$$
      where $M(s_i, s_j)$ is a Matern covariance function and $\phi(x) = c_1 \phi_1(x) + \dots + c_p \phi_p(x)$. 
      
      Here, $\phi_r(x)$ is the $r$th basis function (of a total of $p$) evaluated at $x$, and $c_r,\ r=1, \dots p$ are nonstationary variance parameters
      
      - [Why are these nonstationary?]{.orange}
      
      
      ## Paciorek and Schervish 2006: examples
      
      ```{r echo=FALSE, fig.asp=0.5}
      
      mycppcode <- "
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

arma::mat sigmafy(const arma::rowvec& coordinate, const arma::vec& param){
  return param(0) * arma::diagmat(param(1)*abs(coordinate) + param(2));
}

//[[Rcpp::export]]
arma::mat mycovariance(const arma::mat& cx, const arma::mat& cy, const arma::vec& param){
  arma::mat result = arma::zeros(cx.n_rows, cy.n_rows);
  for(unsigned int i=0; i<cx.n_rows; i++){
    for(unsigned int j=0; j<cy.n_rows; j++){
    
      arma::mat Sigmai = sigmafy(cx.row(i), param);
      arma::mat Sigmaj = sigmafy(cy.row(j), param);
      arma::mat SSinv = arma::inv_sympd(0.5*Sigmai + 0.5*Sigmaj);
      
      double deti = pow(arma::det(Sigmai), -0.25);
      double detj = pow(arma::det(Sigmaj), -0.25);
      double detij = pow(arma::det(0.5*Sigmai + 0.5*Sigmaj), -0.5);
      
      arma::rowvec h = cx.row(i) - cy.row(j);
      double inner = arma::conv_to<double>::from( sqrt(h * SSinv * arma::trans(h)) );
      result(i, j) = deti * detj * detij * exp(-inner);
    }
  }
  return result;
}
"
Rcpp::sourceCpp(code=mycppcode)

coords <- expand.grid(xgrid <- seq(-2, 2, length.out=114), ygrid <- seq(-1,1,length.out=57))
colnames(coords) <- c("xcoord", "ycoord")
cmat <- coords %>% as.matrix()
C <- mycovariance(cmat, cmat, c(1, 0, 1))
L <- t(chol(C))
set.seed(696)
Y <- L %*% rnorm(nrow(coords))

df <- data.frame(coords, Y)

(data_plot <- ggplot(df, aes(xcoord, ycoord, fill=Y)) + 
    geom_raster() +
    scale_fill_scico(palette="bamO") +
    theme_minimal())
```


## Paciorek and Schervish 2006: simulated data examples

```{r echo=FALSE, fig.asp=0.5}
C <- mycovariance(cmat, cmat, c(1, 2, 1))
L <- t(chol(C))
set.seed(696)
Y <- L %*% rnorm(nrow(coords))

df <- data.frame(coords, Y)

(data_plot <- ggplot(df, aes(xcoord, ycoord, fill=Y)) + 
    geom_raster() +
    scale_fill_scico(palette="bamO") +
    theme_minimal())
```

## Paciorek and Schervish 2006: simulated data examples

```{r echo=FALSE, fig.asp=0.5}
C <- mycovariance(cmat, cmat, c(1, 5, 1))
L <- t(chol(C))
set.seed(696)
Y <- L %*% rnorm(nrow(coords))

df <- data.frame(coords, Y)

(data_plot <- ggplot(df, aes(xcoord, ycoord, fill=Y)) + 
    geom_raster() +
    scale_fill_scico(palette="bamO") +
    theme_minimal())
```


## [*896/Advanced topics*] Extending GPs to multivariate data

-  We have considered GPs for univariate data
-  Univariate data: one random variable for each spatial location
-  Multivariate data: one random *vector* for each spatial location
-  Nothing in our discussion changes when we want to define a multivariate GP
-  **Except:** we need to define a *cross*-covariance function
-  Cross-covariance function: not only defines covariance in space, but also across variables
-  $C(\ell_i, \ell_j)$ is a *matrix valued* function. The $(r,s)$ element of $C(\ell_i, \ell_j)$ is the covariance between the $r$ variable measured at location $\ell_i$ and the $s$ variable measured at location $\ell_j$
  -  See book **Chapter 9.2**.
-  More on this at the end of the semester.


## First steps in modeling data via GPs

-  Probability models for point-referenced spatial data interpret observations as realizations of a set of spatially-indexed random variables
-  At each observed spatial location, we "see" 1 realization of a random variable
-  All random variables are dependent on each other
-  We *assume* that dependence generally decays with distance
-  We *assume* that a multivariate Gaussian distribution is appropriate for our collection of random variables
-  We *assume* a covariance model with a small number of parameters
-  We aim to **learn** about the covariance model using the observed data
-  Everything is in place to apply our simple GP model to point-referenced data
