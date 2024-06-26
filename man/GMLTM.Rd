% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GMLTM.R
\name{GMLTM}
\alias{GMLTM}
\title{The Generalized Multidimensional Latent Trait Model}
\usage{
GMLTM(data, Q, components, iters = 2000, iter_warmup = 1000, chains = 2, seed = 1606)
}
\arguments{
\item{data}{$$n \times p$$ data.frame or data matrix with the individuals in rows and items in columns.}

\item{Q}{$$p \times q$$ Q matrix.}

\item{components}{List of $$d$$ elements relating each components to a vector of rules.}

\item{iters}{Number of samples from the posterior distribution.}

\item{chains}{Number of Markov chains.}

\item{iter_warmup}{Number of samples to discard then initializing a Markov chain.}

\item{interval}{Probability associated with the credible intervals.}

\item{seed}{Random seed.}

\item{parallel_chains}{Number of chains to run in parallel cores.}

\item{threads_per_chain}{Number of cores to run within a chain.}

\item{...}{Additional arguments to pass to the \code{sample} function from \code{cmdstanr}.}
}
\value{
\code{cmdstan} object:
}
\description{
Estimate the parameters of the GMLTM.
}
\details{
\code{GMLTM} estimates a ...
}
\references{
Ramirez E., Jiménez M., Franco V., Alvarado J. (2024). Delving into... preprint.
}
