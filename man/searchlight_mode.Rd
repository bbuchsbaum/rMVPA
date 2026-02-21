\name{searchlight_mode}
\alias{searchlight_mode}
\title{Get Or Set Searchlight Mode}
\usage{
searchlight_mode(mode = NULL)
}
\arguments{
\item{mode}{Optional mode value. One of \code{"fast"} or \code{"legacy"}.
If \code{NULL}, returns current effective mode.}
}
\value{
If \code{mode} is \code{NULL}, returns current mode as character.
Otherwise, sets the mode and returns it invisibly.
}
\description{
Simplified high-level runtime switch for searchlight performance behavior.
\code{"fast"} is the default optimized mode and \code{"legacy"} restores
conservative compatibility behavior.
}
