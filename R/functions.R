#trim bulk and single-cell data to contain the same genes
trimData<-function(Signature,bulkData){
	Genes<-intersect(rownames(Signature),names(bulkData))
  	B<-bulkData[Genes]
  	S<-Signature[Genes,,drop=FALSE]
  	return(list("sig"=S,"bulk"=B))
}

#' Solver using OLS, constrained such that cell type numbers>0
#'
#' @param S Signature matrix
#' @param B One column of the bulk data matrix, so one sample
#' @param verbose Whether to produce an output on the console.
#'
#' @return A vector with the cell type proportions for the sample
#'
solveOLS <- function(S, B, verbose = FALSE) {
  solution <- solveOLSInternal(S, B, verbose)
  return(solution / sum(solution))
}

#' Solver using OLS, constrained such that cell type numbers>0
#'
#' Returns absolute numbers, not proportions
#'
#' @param S Signature matrix
#' @param B One column of the bulk data matrix, so one sample
#' @param verbose Whether to produce an output on the console.
#'
#' @return A vector with the cell type numbers for the sample
#'
solveOLSInternal <- function(S, B, verbose = FALSE) {
  D <- t(S) %*% S
  d <- t(S) %*% B
  A <- cbind(diag(dim(S)[2]))
  bzero <- c(rep(0, dim(S)[2]))
  sc <- norm(D, "2")

  solution <- tryCatch(
    {
      quadprog::solve.QP(D / sc, d / sc, A, bzero)$solution
    },
    error = function(cond) {
      message(
        "Error solving the quadratic programming problem. Your dataset might be too small. Run ",
        "with verbose=TRUE to get the full error message."
      )
      if (verbose) {
        stop(cond)
      } else {
        stop()
      }
    },
    warning = function(cond) {
      warning(cond)
    }
  )

  names(solution) <- colnames(S)
  return(solution)
}


#' Solver using WLS with weights dampened by a certain dampening constant
#'
#' @param S Signature matrix
#' @param B One column of the bulk data matrix, so one sample
#' @param verbose Whether to produce an output on the console.
#'
#' @return A vector with the cell type numbers for the sample
#'
solveDampenedWLS <- function(S, B, verbose = FALSE) {
  # first solve OLS, use this solution to find a starting point for the weights
  solution <- solveOLSInternal(S, B, verbose)
  # now use dampened WLS, iterate weights until convergence
  iterations <- 0
  changes <- c()
  # find dampening constant for weights using cross-validation
  j <- findDampeningConstant(S, B, solution, verbose)
  change <- 1
  while (change > .01 & iterations < 1000) {
    newsolution <- solveDampenedWLSj(S, B, solution, j, verbose)
    # decrease step size for convergence
    solutionAverage <-
      rowMeans(cbind(newsolution, matrix(
        solution,
        nrow = length(solution), ncol = 4
      )))
    change <- norm(as.matrix(solutionAverage - solution))
    solution <- solutionAverage
    iterations <- iterations + 1
    changes <- c(changes, change)
  }
  if (verbose) {
    print(round(solution / sum(solution), 5))
  }
  return(solution / sum(solution))
}


#' Solve WLS given a dampening constant
#'
#' @param S Signature matrix
#' @param B One column of the bulk data matrix, so one sample
#' @param goldStandard The average of all the solutions so far
#' @param j The dampening constant
#' @param verbose Whether to produce an output on the console
#'
#' @return A vector with the cell type numbers for the sample
#'
solveDampenedWLSj <-
  function(S, B, goldStandard, j, verbose = FALSE) {
    multiplier <- 1 * 2^(j - 1)
    sol <- goldStandard
    ws <- as.vector((1 / (S %*% sol))^2)
    wsScaled <- ws / min(ws)
    wsDampened <- wsScaled
    wsDampened[which(wsScaled > multiplier)] <- multiplier
    W <- diag(wsDampened)
    D <- t(S) %*% W %*% S
    d <- t(S) %*% W %*% B
    A <- cbind(diag(dim(S)[2]))
    bzero <- c(rep(0, dim(S)[2]))
    sc <- norm(D, "2")

    solution <- tryCatch(
      {
        quadprog::solve.QP(D / sc, d / sc, A, bzero)$solution
      },
      error = function(cond) {
        message(
          "Error solving the quadratic programming problem. Your dataset might be too small. Run ",
          "with verbose=TRUE to get the full error message."
        )
        if (verbose) {
          stop(cond)
        } else {
          stop()
        }
      },
      warning = function(cond) {
        warning(cond)
      }
    )

    names(solution) <- colnames(S)
    return(solution)
  }

#' Finding a dampening constant for the weights using cross-validation
#'
#' @param S Signature matrix
#' @param B One column of the bulk data matrix, so one sample
#' @param goldStandard The solution found with OLS
#' @param verbose Whether to produce an output on the console
#'
#' @return The dampening constant (integer)
findDampeningConstant <- function(S, B, goldStandard, verbose = FALSE) {
  solutionsSd <- NULL
  # goldStandard is used to define the weights
  sol <- goldStandard
  ws <- as.vector((1 / (S %*% sol))^2)
  wsScaled <- ws / min(ws)
  wsScaledMinusInf <- wsScaled
  # ignore infinite weights
  if (max(wsScaled) == "Inf") {
    wsScaledMinusInf <- wsScaled[-which(wsScaled == "Inf")]
  }
  # try multiple values of the dampening constant (multiplier)
  # for each, calculate the variance of the dampened weighted solution for a subset of genes
  for (j in 1:ceiling(log2(max(wsScaledMinusInf)))) {
    multiplier <- 1 * 2^(j - 1)
    wsDampened <- wsScaled
    wsDampened[which(wsScaled > multiplier)] <- multiplier
    solutions <- NULL
    seeds <- c(1:100)

    for (i in 1:100) {
      set.seed(seeds[i]) # make nondeterministic
      subset <-
        sample(length(ws), size = length(ws) * 0.5) # randomly select half of gene set
      # solve dampened weighted least squares for subset
      fit <- stats::lm(B[subset] ~ -1 + S[subset, , drop = FALSE], weights = wsDampened[subset])
      sol <- fit$coef * sum(goldStandard) / sum(fit$coef)
      solutions <- cbind(solutions, sol)
    }

    solutionsSd <-
      cbind(solutionsSd, apply(solutions, 1, stats::sd))
  }
  # choose dampening constant that results in least cross-validation variance
  j <- which.min(colMeans(solutionsSd^2))
  return(j)
}


#' Solver using SVR
#'
#' @param S Signature matrix
#' @param B One column of the bulk data matrix, so one sample
#' @param verbose Whether to produce an output on the console
#'
#' @return A vector with the cell type proportions for the sample
#'
solveSVR <- function(S, B, verbose = FALSE) {
  # scaling
  ub <- max(c(as.vector(S), B)) # upper bound
  lb <- min(c(as.vector(S), B)) # lower bound
  Bs <- ((B - lb) / ub) * 2 - 1
  Ss <- ((S - lb) / ub) * 2 - 1

  # perform SVR
  model <-
    e1071::svm(
      Ss,
      Bs,
      nu = 0.5,
      scale = TRUE,
      type = "nu-regression",
      kernel = "linear",
      cost = 1
    )
  coef <- t(model$coefs) %*% model$SV
  coef[which(coef < 0)] <- 0
  coef <- as.vector(coef)
  names(coef) <- colnames(S)
  if (verbose) {
    print(round(coef / sum(coef), 5))
  }
  return(coef / sum(coef))
}


#' Performing DE analysis using Seurat
#'
#' When path = NULL, the generated files in the processes will not be saved and output
#'
#' @param scdata The single cell data matrix
#' @param id A Vector of the cell type annotations
#' @param path OPTIONAL path for saving generated files
#' @param verbose Whether to produce an output on the console
#'
#' @return List with the differentially expressed genes for each cell type
#'
DEAnalysis <- function(scdata, id, path = NULL, verbose = FALSE) {
  uniqueIds <- unique(id)
  list.DE.group <- as.list(rep(0, length(uniqueIds)))

  expr_obj <-
    Seurat::CreateSeuratObject(counts = as.data.frame(scdata), project = "DE")
  expr_obj2 <- Seurat::SetIdent(expr_obj, value = as.vector(id))
  if (verbose) {
    print("Calculating differentially expressed genes:")
  }
  for (i in uniqueIds) {
    # check if there are more than 3 cells in the current id group
    if(length(which(expr_obj2@active.ident == i)) < 3){
      message(paste0('There are fewer than 3 cells in cell type ', i,
                     '. Please remove this cell type, otherwise no differential expression analysis can be performed using Seurat.'))
      stop()
    }
    de_group <-
      Seurat::FindMarkers(
        object = expr_obj2,
        ident.1 = i,
        ident.2 = NULL,
        only.pos = TRUE,
        test.use = "bimod"
      )

    index <- which(uniqueIds == i)
    list.DE.group[[index]] <- de_group

    if (!is.null(path)) {
      save(de_group, file = paste(path, "/de_", i, ".RData", sep = ""))
    }
  }
  return(list.DE.group)
}


#' Building the signature matrix using Seurat
#'
#' When path = NULL, the generated files in the processes will not be saved and output.
#'
#' @param scdata The single cell data matrix
#' @param id A Vector of the cell type annotations
#' @param path OPTIONAL path for saving generated files
#' @param verbose Whether to produce an output on the console
#' @param diff.cutoff The FC cutoff
#' @param pval.cutoff The pValue cutoff
#'
#' @return The computed signature matrix
#'
buildSignatureMatrixUsingSeurat <- function(scdata,
                                            id,
                                            path = NULL,
                                            verbose = FALSE,
                                            diff.cutoff = 0.5,
                                            pval.cutoff = 0.01) {
  # perform differential expression analysis
  list.DE.groups <- DEAnalysis(scdata, id, path, verbose)

  numberofGenes <- c()
  uniqueIds <- unique(id)
  for (i in uniqueIds) {
    if (file.exists(paste(path, "/de_", i, ".RData", sep = ""))) {
      load(file = paste(path, "/de_", i, ".RData", sep = ""))
    } else {
      index <- which(uniqueIds == i)
      de_group <- list.DE.groups[[index]]
    }

    DEGenes <-
      rownames(de_group)[intersect(
        which(de_group$p_val_adj < pval.cutoff),
        which(de_group$avg_log2FC > diff.cutoff)
      )]
    nonMir <- grep("MIR|Mir", DEGenes, invert = TRUE)
    assign(
      paste("cluster_lrTest.table.", i, sep = ""),
      de_group[which(rownames(de_group) %in% DEGenes[nonMir]), ]
    )
    numberofGenes <- c(numberofGenes, length(DEGenes[nonMir]))

    if(numberofGenes == 0){
      message('Did not find any significant marker genes for cell type ', i,'. Please adjust cutoff parameters.')
      stop()
    }
  }

  # need to reduce number of genes
  # for each subset, order significant genes by decreasing fold change,
  # choose between 50 and 200 genes
  # choose matrix with lowest condition number
  conditionNumbers <- c()
  for (g in 50:200) {
    Genes <- c()
    j <- 1
    for (i in uniqueIds) {
      if (numberofGenes[j] > 0) {
        temp <- paste("cluster_lrTest.table.", i, sep = "")
        temp <- get(temp)
        temp <- temp[order(temp$p_val_adj, decreasing = TRUE), ]
        Genes <-
          c(Genes, (rownames(temp)[1:min(g, numberofGenes[j])]))
      }
      j <- j + 1
    }
    Genes <- unique(Genes)
    # make signature matrix
    ExprSubset <- scdata[Genes, , drop = FALSE]
    Sig <- NULL
    for (i in uniqueIds) {
      Sig <-
        cbind(Sig, (apply(ExprSubset, 1, function(y) {
          mean(y[which(id == i)])
        })))
    }
    colnames(Sig) <- uniqueIds
    conditionNumbers <- c(conditionNumbers, kappa(Sig))
  }

  # g is optimal gene number
  g <- which.min(conditionNumbers) + min(49, numberofGenes - 1)
  Genes <- c()
  j <- 1
  for (i in uniqueIds) {
    if (numberofGenes[j] > 0) {
      temp <- paste("cluster_lrTest.table.", i, sep = "")
      temp <- get(temp)
      temp <- temp[order(temp$p_val_adj, decreasing = TRUE), ]
      Genes <-
        c(Genes, (rownames(temp)[1:min(g, numberofGenes[j])]))
    }
    j <- j + 1
  }
  Genes <- unique(Genes)
  ExprSubset <- scdata[Genes, , drop = FALSE]
  Sig <- NULL
  for (i in uniqueIds) {
    Sig <-
      cbind(Sig, (apply(ExprSubset, 1, function(y) {
        mean(y[which(id == i)])
      })))
  }

  colnames(Sig) <- uniqueIds
  rownames(Sig) <- Genes

  if (!is.null(path)) {
    save(Sig, file = paste(path, "/Sig.RData", sep = ""))
  }

  return(Sig)
}


# Functions needed for DE

Mean.in.log2space <- function(x, pseudo.count) {
  return(log2(mean(2^(x) - pseudo.count) + pseudo.count))
}

stat.log2 <- function(data.m, group.v, pseudo.count) {
  # data.m=data.used.log2
  log2.mean.r <-
    stats::aggregate(t(data.m), list(as.character(group.v)), function(x) {
      Mean.in.log2space(x, pseudo.count)
    })
  log2.mean.r <- t(log2.mean.r)
  colnames(log2.mean.r) <-
    paste("mean.group", log2.mean.r[1, ], sep = "")
  log2.mean.r <- log2.mean.r[-1, ]
  log2.mean.r <- as.data.frame(log2.mean.r)
  log2.mean.r <- varhandle::unfactor(log2.mean.r) # from varhandle
  log2.mean.r[, 1] <- as.numeric(log2.mean.r[, 1])
  log2.mean.r[, 2] <- as.numeric(log2.mean.r[, 2])
  log2_foldchange <- log2.mean.r$mean.group1 - log2.mean.r$mean.group0
  results <- data.frame(cbind(
    log2.mean.r$mean.group0,
    log2.mean.r$mean.group1,
    log2_foldchange
  ))
  colnames(results) <- c("log2.mean.group0", "log2.mean.group1", "log2_fc")
  rownames(results) <- rownames(log2.mean.r)
  return(results)
}

v.auc <- function(data.v, group.v) {
  prediction.use <- ROCR::prediction(data.v, group.v, 0:1)
  perf.use <- ROCR::performance(prediction.use, "auc")
  auc.use <- round(perf.use@y.values[[1]], 3)
  return(auc.use)
}

m.auc <- function(data.m, group.v) {
  AUC <- unlist(mclapply(seq(1, nrow(data.m)), function(i) {
    x <- unlist(data.m[i,])
    v.auc(x, group.v)
  }))
  AUC[is.na(AUC)] <- 0.5
  return(AUC)
}


#' Performing DE analysis using mast
#'
#' When path = NULL, the generated files in the processes will not be saved and output.
#'
#' @param scdata The single cell data matrix
#' @param id A Vector of the cell type annotations
#' @param path OPTIONAL path for saving generated files
#' @param verbose Whether to produce an output on the console.
#'
#' @return A list with the cell types and their differentially expressed genes
#'
DEAnalysisMAST <- function(scdata, id, path, verbose = FALSE) {
  uniqueIds <- unique(id)
  list_lrTest.table <- as.list(rep(0, length(uniqueIds)))

  pseudo.count <- 0.1
  data.used.log2 <- log2(scdata + pseudo.count)
  colnames(data.used.log2) <- make.unique(colnames(data.used.log2))
  diff.cutoff <- 0.5
  for (i in uniqueIds) {
    cells.symbol.list2 <- colnames(data.used.log2)[which(id == i)]
    cells.coord.list2 <- match(cells.symbol.list2, colnames(data.used.log2))
    cells.symbol.list1 <- colnames(data.used.log2)[which(id != i)]
    cells.coord.list1 <- match(cells.symbol.list1, colnames(data.used.log2))
    data.used.log2.ordered <-
      cbind(data.used.log2[, cells.coord.list1], data.used.log2[, cells.coord.list2])
    group.v <-
      c(rep(0, length(cells.coord.list1)), rep(1, length(cells.coord.list2)))
    # ouput
    log2.stat.result <-
      stat.log2(data.used.log2.ordered, group.v, pseudo.count)
    Auc <- m.auc(data.used.log2.ordered, group.v)
    bigtable <- data.frame(cbind(log2.stat.result, Auc))

    de <- bigtable[bigtable$log2_fc > diff.cutoff, ]
    if (verbose) {
      dim(de)
    }
    if (dim(de)[1] > 1) {
      data.1 <- data.used.log2[, cells.coord.list1, drop = FALSE]
      data.2 <- data.used.log2[, cells.coord.list2, drop = FALSE]
      genes.list <- rownames(de)
      log2fold_change <- cbind(genes.list, de$log2_fc)
      colnames(log2fold_change) <- c("gene.name", "log2fold_change")
      counts <- as.data.frame(cbind(data.1[genes.list, ], data.2[genes.list, ]))
      groups <- c(
        rep("Cluster_Other", length(cells.coord.list1)),
        rep(i, length(cells.coord.list2))
      )
      groups <- as.character(groups)
      data_for_MIST <-
        verbose_wrapper(verbose)(as.data.frame(cbind(
          rep(rownames(counts), dim(counts)[2]),
          reshape::melt(counts),
          rep(groups, each = dim(counts)[1]),
          rep(1, dim(counts)[1] * dim(counts)[2])
        )))
      colnames(data_for_MIST) <- c(
        "gene",
        "Subject.ID",
        "Et",
        "Population",
        "Number.of.Cells"
      )
      vbeta <- data_for_MIST
      vbeta.fa <-
        verbose_wrapper(verbose)(
          MAST::FromFlatDF(
            vbeta,
            idvars = c("Subject.ID"),
            primerid = "gene",
            measurement = "Et",
            ncells = "Number.of.Cells",
            geneid = "gene",
            cellvars = c("Number.of.Cells", "Population"),
            phenovars = c("Population"),
            id = "vbeta all"
          )
        )
      vbeta.1 <- subset(vbeta.fa, Number.of.Cells == 1)
      # .3 mast
      # utils::head(SummarizedExperiment::colData(vbeta.1)) ??
      zlm.output <-
        verbose_wrapper(verbose)(MAST::zlm(
          ~Population,
          vbeta.1,
          method = "bayesglm",
          ebayes = TRUE
        ))
      if (verbose) {
        methods::show(zlm.output)
        coefAndCI <- summary(zlm.output, logFC = TRUE)
        print(coefAndCI)
      }
      zlm.lr <-
        verbose_wrapper(verbose)(MAST::lrTest(zlm.output, "Population"))
      zlm.lr_pvalue <- reshape::melt(zlm.lr[, , "Pr(>Chisq)"])
      zlm.lr_pvalue <-
        zlm.lr_pvalue[which(zlm.lr_pvalue$test.type == "hurdle"), ]



      lrTest.table <-
        merge(zlm.lr_pvalue, de, by.x = "primerid", by.y = "row.names")
      colnames(lrTest.table) <-
        c(
          "gene",
          "test.type",
          "p_value",
          paste("log2.mean.", "Cluster_Other", sep = ""),
          paste("log2.mean.", i, sep = ""),
          "log2fold_change",
          "Auc"
        )
      cluster_lrTest.table <-
        lrTest.table[rev(order(lrTest.table$Auc)), ]

      # . 4 save results

      index <- which(uniqueIds == i)
      list_lrTest.table[[index]] <- cluster_lrTest.table

      if (!is.null(path)) {
        utils::write.csv(cluster_lrTest.table,
                         file = paste(path, "/", i, "_lr_test.csv", sep = "")
        )
        save(cluster_lrTest.table,
             file = paste(path, "/", i, "_mist.RData", sep = "")
        )
      }
    }
  }

  return(list_lrTest.table)
}


#' #' Building the signature matrix using mast
#'
#' When path = NULL, the generated files in the processes will not be saved and output.
#'
#' @param scdata The single cell data matrix
#' @param id A Vector of the cell type annotations
#' @param path OPTIONAL path for saving generated files
#' @param verbose Whether to produce an output on the console.
#' @param ncores How many cores to use for DGE analysis
#' @param diff.cutoff The FC cutoff
#' @param pval.cutoff The pValue cutoff
#'
#' @return The computed signature matrix
#' 
#' @import parallel
#'
buildSignatureMatrixMAST <- function(scdata,
                                     id,
                                     path = NULL,
                                     verbose = FALSE,
                                     ncores = 1,
                                     diff.cutoff = 0.5,
                                     pval.cutoff = 0.01) {
  # number of cores for:
  # m.auc
  # MAST functions
  options(mc.cores=ncores)

  # compute differentially expressed genes for each cell type
  list.cluster.table <-
    DEAnalysisMAST(scdata, id, path, verbose = verbose)

  # for each cell type, choose genes in which FDR adjusted p-value is less than 0.01 and the
  # estimated fold-change is greater than 0.5
  numberofGenes <- c()
  uniqueIds <- unique(id)
  for (i in uniqueIds) {
    if (file.exists(paste(path, "/", i, "_MIST.RData", sep = ""))) {
      load(file = paste(path, "/", i, "_MIST.RData", sep = ""))
    } else {
      index <- which(uniqueIds == i)
      cluster_lrTest.table <- list.cluster.table[[index]]
    }
    pvalue_adjusted <-
      stats::p.adjust(
        cluster_lrTest.table[, 3],
        method = "fdr",
        n = length(cluster_lrTest.table[, 3])
      )
    cluster_lrTest.table <-
      cbind(cluster_lrTest.table, pvalue_adjusted)
    DEGenes <-
      cluster_lrTest.table$gene[intersect(
        which(pvalue_adjusted < pval.cutoff),
        which(cluster_lrTest.table$log2fold_change > diff.cutoff)
      )]

    # because Mir gene is usually not accurate
    nonMir <- grep("MIR|Mir", DEGenes, invert = TRUE)
    assign(
      paste("cluster_lrTest.table.", i, sep = ""),
      cluster_lrTest.table[which(cluster_lrTest.table$gene %in% DEGenes[nonMir]), ]
    )
    numberofGenes <- c(numberofGenes, length(DEGenes[nonMir]))
  }


  # need to reduce number of genes
  # for each subset, order significant genes by decreasing fold change,
  # choose between 50 and 200 genes
  # for each, iterate and choose matrix with lowest condition number
  conditionNumbers <- c()
  for (g in 50:200) {
    Genes <- c()
    j <- 1
    for (i in uniqueIds) {
      if (numberofGenes[j] > 0) {
        temp <- paste("cluster_lrTest.table.", i, sep = "")
        temp <- get(temp)
        temp <-
          temp[order(temp$log2fold_change, decreasing = TRUE), ]
        Genes <-
          c(Genes, varhandle::unfactor(temp$gene[1:min(g, numberofGenes[j])]))
      }
      j <- j + 1
    }
    Genes <- unique(Genes)
    # make signature matrix
    ExprSubset <- scdata[Genes, , drop = FALSE]
    Sig <- NULL
    for (i in uniqueIds) {
      Sig <-
        cbind(Sig, (apply(ExprSubset, 1, function(y) {
          mean(y[which(id == i)])
        })))
    }
    colnames(Sig) <- uniqueIds
    conditionNumbers <- c(conditionNumbers, kappa(Sig))
  }
  # g is optimal gene number
  g <- which.min(conditionNumbers) + min(49, numberofGenes - 1)
  Genes <- c()
  j <- 1
  for (i in uniqueIds) {
    if (numberofGenes[j] > 0) {
      temp <- paste("cluster_lrTest.table.", i, sep = "")
      temp <- get(temp)
      temp <- temp[order(temp$log2fold_change, decreasing = TRUE), ]
      Genes <-
        c(Genes, varhandle::unfactor(temp$gene[1:min(g, numberofGenes[j])]))
    }
    j <- j + 1
  }
  Genes <- unique(Genes)
  ExprSubset <- scdata[Genes, , drop = FALSE]
  Sig <- NULL
  for (i in uniqueIds) {
    Sig <-
      cbind(Sig, (apply(ExprSubset, 1, function(y) {
        mean(y[which(id == i)])
      })))
  }

  colnames(Sig) <- uniqueIds
  rownames(Sig) <- Genes

  if (!is.null(path)) {
    save(Sig, file = paste(path, "/Sig.RData", sep = ""))
  }


  return(Sig)
}

#' A wrapper function whether to suppress messages
#'
#' @param verbose Whether to produce an output on the console.
#'
#' @return A function which will suppress messages or not, depending on the verbose parameter
#'
verbose_wrapper <- function(verbose) {
  return(function(method) {
    if (!verbose) {
      suppressMessages(method)
    } else {
      method
    }
  })
}
