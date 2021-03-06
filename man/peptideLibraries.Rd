% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/library.R
\name{peptideLibraries}
\alias{peptideLibraries}
\alias{getLibraryPath}
\alias{setLibraryPath}
\alias{getLibrary}
\alias{makeLibrary}
\title{Peptide libraries}
\usage{
getLibraryPath()

setLibraryPath(path)

getLibrary(name)

makeLibrary(library, name)
}
\arguments{
\item{path}{path to a folder containing a .rds files for each library}

\item{name}{name of the library}

\item{library}{a \code{matrix}, \code{data.frame}, or
\linkS4class{DataFrame} with the peptide information for the
specified library.}
}
\value{
\code{getLibraryPath} returns the path to the directory containing
.rds files for each library. \code{getLibrary} returns a
\linkS4class{DataFrame} corresponding to the peptide information for the
specified library.
}
\description{
PhIP-Seq experiments often use identical peptide libraries
different cohorts. These functions enable the user to conveniently reuse
tidied libraries.
}
\details{
Each library is stored as a \linkS4class{DataFrame} in .rds file.
By default the libraries are stored in the \code{libraries} package folder.
\code{setLibraryPath} allows the user to change the location of this library,
if desired.

While libraries can be directly saved in the specified format by the user,
it is highly recommended to use the \code{makeLibrary} function to
save a new library.
}
\section{Functions}{
\itemize{
\item \code{getLibraryPath}: return the path to a folder containing the
libraries.

\item \code{setLibraryPath}: set the path to a folder containing the
libraries.

\item \code{getLibrary}: return a \linkS4class{DataFrame} with the
peptide information corresponding to the library.

\item \code{makeLibrary}: create and store a \linkS4class{DataFrame} with
the specified peptide information.
}}

\examples{
## Get and set path to libraries folder
getLibraryPath()
\dontrun{setLibraryPath("examplepath/")}

## Create a new library
pep_meta <- data.frame(species = c(rep("human immunodeficiency virus", 3),
     rep("Epstein-Barr virus", 2)))
makeLibrary(pep_meta, "new_library")

## Use new library
counts_dat <- matrix(1:10, nrow = 5)
phip_obj <- PhIPData(counts = counts_dat,
     peptideInfo = getLibrary("new_library"))

## Delete created library
file.remove(paste0(getLibraryPath(), "/new_library.rds"))

}
