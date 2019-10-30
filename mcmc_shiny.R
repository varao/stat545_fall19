# MCMC Shiny demo for STAT545
# Vinayak Rao

library(ggplot2)
library(shiny)

a_par <- .3
b_par <- 3

rosenbrock <- function(x,y,a=a_par,b=b_par) {
  val <- (a-x)^2 + b*(y/1-x*x)^2
  return(val)
}

rosenbrock_grad <- function(x,y,a=a_par,b=b_par) {
  val <- c(-2*(a-x) - 4*b*x*(y-x*x),
           2*b*(y-x*x)
           )
  return(val)
}
# Modified Neal's code to suit my purposes: Vinayak Rao
# SIMPLE IMPLEMENTATION OF HAMILTONIAN MONTE CARLO.
#
# Radford M. Neal, 2010.
#
# This program appears in Figure 2 of "MCMC using Hamiltonian dynamics",
# to appear in the Handbook of Markov Chain Monte Carlo.
#
# The arguments to the HMC function are as follows:
#
#   U          A function to evaluate minus the log of the density of the
#              distribution to be sampled, plus any constant - ie, the
#              "potential energy".
#
#   grad_U     A function to evaluate the gradient of U.
#
#   epsilon    The stepsize to use for the leapfrog steps.
#
#   L          The number of leapfrog steps to do to propose a new state.
#
#   current_q  The current state (position variables only).
#
# Momentum variables are sampled from independent standard normal
# distributions within this function.  The value return is the vector
# of new position variables (equal to current_q if the endpoint of the
# trajectory was rejected).
#
# This function was written for illustrative purposes.  More elaborate
# implementations of this basic HMC method and various variants of HMC
# are available from my web page, http://www.cs.utoronto.ca/~radford/


HMC = function (U, grad_U, L, current_q) {
  epsilon = scl
  L       = run
  M       = diag(c(1,1))
  M       = matrix(c(1,.8,.8,1),nrow=2)
  #M       = matrix(c(1,-.8,-.8,1),nrow=2)
  cM      = chol(M)
  iM      = chol2inv(cM)

  q = current_q
  p = cM %*% rnorm(length(q),0,1)  # independent standard normal variates
  current_p = p
  q_vec = matrix(0,L+1,length(q))
  q_vec[1,] = q

  # Make a half step for momentum at the beginning
  p = p - epsilon * grad_U(q[1],q[2]) / 2

  # Alternate full steps for position and momentum
  for (i in 1:L) {
    # Make a full step for the position
    q = q + epsilon * iM %*% p
    q_vec[i+1,] = q
    # Make a full step for the momentum, except at end of trajectory
    if (i!=L) p = p - epsilon * grad_U(q[1],q[2])
  }
  # Make a half step for momentum at the end.
  p = p - epsilon * grad_U(q[1],q[2]) / 2

  # Negate momentum at end of trajectory to make the proposal symmetric
  p = -p

  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(current_q[1],current_q[2])
  current_K = (t(current_p)%*%iM%*%current_p) / 2
  proposed_U = U(q[1],q[2])
  proposed_K = t(p) %*% iM %*% p / 2

  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position

  if(!any(is.nan(q)) &&
    log(runif(1)) < (current_U-proposed_U+current_K-proposed_K)) {
    return (list(rslt=q,path=data.frame(x=q_vec[,1],y=q_vec[,2])))  # accept
  } else {
    return (list(rslt=current_q,path=data.frame(x=q_vec[,1],y=q_vec[,2])))  # accept
  }
}

x <- seq(-2,2.5,.05)
y <- seq(-2,5.5,.05)
df <- data.frame(x=rep(x,length(y)), y = rep(y,each=length(x)))
df$z <- -rosenbrock(df$x,df$y)
rosen <- ggplot(df) + geom_tile(aes(x=x,y=y,fill=exp(z))) +
    coord_cartesian(xlim=c(-2,2.5),  ylim=c(-3,6))

ii <- 1
scl <- 1
run <- 10
pts <- data.frame(x=rep(0,1000), y=rep(0,1000))
mcmc_type  <- 1

runApp(list(
# ui = pageWithSidebar(
  ui = fluidPage(

    headerPanel("MCMC on the Rosenbrock density"),

    selectInput("mcmc_type", label = h3("MCMC algorithm"),
           choices = list("Metropolis-Hastings" = 1, "Gibbs" = 2, "HMC"=3),
           selected = 3),

    sidebarPanel(
      sliderInput("len_scl",
                  "Step size:",
                  min = .01,
                  max = 2,
                  value = .1),
      sliderInput("len_run",
                  "Number of leapfrog steps:",
                  min = 1,
                  max = 200,
                  value = 20)
    ),
    mainPanel(
      plotOutput("distPlot", click = "plot_click")
    )

#   mainPanel(
#     plotOutput("histPlot")
#   )
  ),
  server =function(input, output, session) {
    autoInvalidate <- reactiveTimer(500, session)

    output$distPlot <- renderPlot({
      autoInvalidate()
      # generate an rnorm distribution and plot it
      ii       <<- ii + 1
      if(mcmc_type == 1) {
        prop     <- pts[ii-1,] + scl * rnorm(2);
        pts[ii,] <<- prop
        if(log(runif(1)) > -rosenbrock(pts[ii,1], pts[ii,2]) + rosenbrock(pts[ii-1,1],pts[ii-1,2]))
          pts[ii,] <<- pts[ii-1,]
        rosen + geom_path(data = pts[1:ii,], aes(x=x,y=y), color='red',size=2, alpha=.4) +
          geom_point(data=prop, aes(x=x,y=y), size=4,color='black')
      } else if(mcmc_type == 2) {
        pts[ii,] <<- pts[ii-1,]
        if(ii%%2) {
          fx     <- function(x) {exp(-rosenbrock(x,pts[ii,2]))}
          z      <- integrate(fx,-Inf,Inf)$value
          fx_cdf <- function(x) {integrate(fx,-Inf,x)$value/z - runif(1)}
          pts[ii,1] <<- uniroot(fx_cdf,lower = -6, upper = 6)$root
        } else {
          pts[ii,2] <<- rnorm(1, pts[ii,1]^2, .5/sqrt(b_par))
        }
        rosen + geom_path(data = pts[1:ii,], aes(x=x,y=y), color='red',size=2, alpha=.4) +
          geom_point(data=pts[ii,], aes(x=x,y=y), size=4,color='black')
      } else {
        prop     <- HMC(rosenbrock, rosenbrock_grad, 200,
                                    c(pts[ii-1,1],pts[ii-1,2]))
        pts[ii,] <<- prop$rslt
        if(pts[ii,1] != pts[ii-1,1]) {
          rosen + geom_path(data = pts[1:ii,], aes(x=x,y=y), color='red',size=2, alpha=.4) +
            geom_path(data=prop$path, aes(x=x,y=y), size=1,color='black',
            arrow=arrow(length=unit(0.3,"cm")))
        } else {
          rosen + geom_path(data = pts[1:ii,], aes(x=x,y=y), color='red',size=2, alpha=.4) +
            geom_path(data=prop$path, aes(x=x,y=y), size=1,color='grey',
            arrow=arrow(length=unit(0.3,"cm")))
        }

      }
    })

    output$histPlot <- renderPlot({
      autoInvalidate()
      ggplot(data = pts[1:ii,]) + geom_histogram(aes(x=x),bins = 30)
    })

   observeEvent(input$plot_click, {
     pts[1,] <<- c(input$plot_click$x, input$plot_click$y)
     ii      <<- 1
   })

    observeEvent(input$len_scl, {
      scl <<- input$len_scl
    })

    observeEvent(input$len_run, {
      run <<- input$len_run
    })

    observeEvent(input$mcmc_type, {
      mcmc_type <<- input$mcmc_type
    })

 }
))

#q_vec <- matrix(0,20,2)
#for(ii in 2:20) {
#  q_vec[ii,] <- HMC(rosenbrock,rosenbrock_grad,.1,5,q_vec[ii-1,])$rslt
#}
