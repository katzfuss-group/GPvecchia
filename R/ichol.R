#' Wrapper for incomplete Cholesky decomposition
#'
#' @param M the matrix to be decomposed
#' @param S sparsity pattern matrix given
#'
#' @return the incomplete Cholesky factor in the sparse format
#' @examples
#' A = matrix(runif(25), ncol = 5)
#' A = t(A) * A + 2 * Matrix::Diagonal(5)
#' S = Matrix::Matrix(c(rep(1, 5), c(0, 1, 1, 0, 0), c(0, 0, 1, 0, 1), c(0, 0, 0, 1, 0), c(0, 0, 0, 0, 1)), ncol = 5, byrow = TRUE)
#' I1 = ichol(A, S)
#' I2 = ichol(A * S)
#'
#' @export
ichol = function(M, S = NULL) {

    #### checking input
    if(dim(M)[1] !=dim(M)[2]) {
        stop("The matrix to be decomposed is not square.")
    }
    
    if (!is.null(S) && sum(abs(dim(M) - dim(S))) > 0) {
        stop("sparsity pattern has different dimensions than the matrix")
    }

    # we need only the upper triangle
    if(!methods::is(M, "CsparseMatrix") || !Matrix::isTriangular(M)) {
        M = methods::as(Matrix::triu(M), "CsparseMatrix")
    }
    
    if(Matrix::isTriangular(M) && attr(Matrix::isTriangular(M), "kind")=="L") {
        M = Matrix::t(M)
    }
    
    if(!is.null(S)) {
        if(!Matrix::isTriangular(S)) {
            S = methods::as(Matrix::triu(S), "CsparseMatrix")
        } else if (attr(Matrix::isTriangular(S), "kind") == "L") {
            S = Matrix::t(S)
        }
    }
    
    if(!methods::is(M, "sparseMatrix") && dim(M)[1] > 100) {
        warning("Decomposing a large dense matrix. This might be slow depending on the sparsity pattern.")
    }
    
    #### Actual decomposition
    if(!is.null(S)){
        p=S@p; i=S@i
        M = M * S
    } else {
        p=M@p; i=M@i
    }
    vals = GPvecchia::ic0(p, i, M@x)
    Msp = Matrix::sparseMatrix(i = i, p = p, x = vals, index1 = FALSE)
    Msp@x = vals
    return(Msp)
}
