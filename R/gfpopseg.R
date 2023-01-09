 #' myGraphFUN
#'
#' @param do.this a string of the form Tj where j = 1, ...., 10 represents the number of breakpoints
#'
#' @return a graph needed to compute j breakpoints using gfpop algorithm
#' @export
#' @importFrom gfpop graph Edge StartEnd
myGraphFUN <- function(do.this){
  switch(do.this,
         T1={myGraph <- graph(Edge("1", "2", type = "std", penalty = 0),
                              StartEnd(start = "1", end = "2"),
                              all.null.edges = TRUE)
         },
         T2={myGraph <- graph(Edge("1", "2", type = "std", penalty = 0),
                              Edge("2", "3", type = "std", penalty = 0),
                              StartEnd(start = "1", end = "3"),
                              all.null.edges = TRUE)
         },
         T3={myGraph <- graph(Edge("1", "2", type = "std", penalty = 0),
                              Edge("2", "3", type = "std", penalty = 0),
                              Edge("3", "4", type = "std", penalty = 0),
                              StartEnd(start = "1", end = "4"),
                              all.null.edges = TRUE)
         },
         T4={myGraph <- graph(Edge("1", "2", type = "std", penalty = 0),
                              Edge("2", "3", type = "std", penalty = 0),
                              Edge("3", "4", type = "std", penalty = 0),
                              Edge("4", "5", type = "std", penalty = 0),
                              StartEnd(start = "1", end = "5"),
                              all.null.edges = TRUE)
         },
         T5={myGraph <- graph(Edge("1", "2", type = "std", penalty = 0),
                              Edge("2", "3", type = "std", penalty = 0),
                              Edge("3", "4", type = "std", penalty = 0),
                              Edge("4", "5", type = "std", penalty = 0),
                              Edge("5", "6", type = "std", penalty = 0),
                              StartEnd(start = "1", end = "6"),
                              all.null.edges = TRUE)
         },
         T6={myGraph <- graph(Edge("1", "2", type = "std", penalty = 0),
                              Edge("2", "3", type = "std", penalty = 0),
                              Edge("3", "4", type = "std", penalty = 0),
                              Edge("4", "5", type = "std", penalty = 0),
                              Edge("5", "6", type = "std", penalty = 0),
                              Edge("6", "7", type = "std", penalty = 0),
                              StartEnd(start = "1", end = "7"),
                              all.null.edges = TRUE)
         },
         T7={myGraph <- graph(Edge("1", "2", type = "std", penalty = 0),
                              Edge("2", "3", type = "std", penalty = 0),
                              Edge("3", "4", type = "std", penalty = 0),
                              Edge("4", "5", type = "std", penalty = 0),
                              Edge("5", "6", type = "std", penalty = 0),
                              Edge("6", "7", type = "std", penalty = 0),
                              Edge("7", "8", type = "std", penalty = 0),
                              StartEnd(start = "1", end = "8"),
                              all.null.edges = TRUE)
         },
         T8={myGraph <- graph(Edge("1", "2", type = "std", penalty = 0),
                              Edge("2", "3", type = "std", penalty = 0),
                              Edge("3", "4", type = "std", penalty = 0),
                              Edge("4", "5", type = "std", penalty = 0),
                              Edge("5", "6", type = "std", penalty = 0),
                              Edge("6", "7", type = "std", penalty = 0),
                              Edge("7", "8", type = "std", penalty = 0),
                              Edge("8", "9", type = "std", penalty = 0),

                              StartEnd(start = "1", end = "9"),
                              all.null.edges = TRUE)
         },
         T9={myGraph <- graph(Edge("1", "2", type = "std", penalty = 0),
                              Edge("2", "3", type = "std", penalty = 0),
                              Edge("3", "4", type = "std", penalty = 0),
                              Edge("4", "5", type = "std", penalty = 0),
                              Edge("5", "6", type = "std", penalty = 0),
                              Edge("6", "7", type = "std", penalty = 0),
                              Edge("7", "8", type = "std", penalty = 0),
                              Edge("8", "9", type = "std", penalty = 0),
                              Edge("9", "10", type = "std", penalty = 0),
                              StartEnd(start = "1", end = "10"),
                              all.null.edges = TRUE)
         },
         T10={myGraph <- graph(Edge("1", "2", type = "std", penalty = 0),
                               Edge("2", "3", type = "std", penalty = 0),
                               Edge("3", "4", type = "std", penalty = 0),
                               Edge("4", "5", type = "std", penalty = 0),
                               Edge("5", "6", type = "std", penalty = 0),
                               Edge("6", "7", type = "std", penalty = 0),
                               Edge("7", "8", type = "std", penalty = 0),
                               Edge("8", "9", type = "std", penalty = 0),
                               Edge("9", "10", type = "std", penalty = 0),
                               Edge("10", "11", type = "std", penalty = 0),

                               StartEnd(start = "1", end = "11"),
                               all.null.edges = TRUE)
         },
         stop("Enter something that switches me!")
  )
  return(myGraph)
}


#' gfpop_gauss
#'
#' @param y a vector of n observations (supposed gaussian and homoscedastic)
#' @param K a parameter specifying that we want K-1 breakpoints
#'
#' @return gfpop results and computed breakpoints
#' @export
gfpop_gauss <- function(y,K){ # gfpop est par défaut homoscedastique !?
  lpi_jk = matPijk(K)
  myGraph <- myGraphFUN(paste("T",K-1,sep = ""))
  gfpop_res <- gfpop(data = y , mygraph = myGraph , type = "mean")
  return(list(res_gfpop = gfpop_res, bp_gfpop = gfpop_res$changepoints[1:(K-1)]))
}
# res_gfpop_gauss =  gfpop_gauss(y,lpi_jk)
# res_gfpop_gauss$bp_gfpop


#' gfpop_poisson
#'
#' @param y a vector of n observations (supposed poissonian and homoscedastic)
#' @param lpi_jk a KxK matrix specifying that we want K-1 breakpoints (to be simplified !!!)
#'
#' @return gfpop results and computed breakpoints
#' @export
#' @importFrom gfpop gfpop
gfpop_poisson <- function(y,lpi_jk){ # gfpop est par défaut homoscedastique !? On peut l'utiliser/adapater pour survreg !?
  K = dim(lpi_jk)[1]
  myGraph <- myGraphFUN(paste("T",K-1,sep = ""))
  gfpop_res <- gfpop(data = y , mygraph = myGraph , type = "poisson")
  return(list(res_gfpop = gfpop_res, bp_gfpop = gfpop_res$changepoints[1:(K-1)]))
}
# res_gfpop_poisson =  gfpop_poisson(y,lpi_jk)
# res_gfpop_poisson$bp_gfpop
