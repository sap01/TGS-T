/*
 * Goal: Calc score of a given Bayesian
 * network w.r.t. a given scoring 
 * function.
 * 
 * This source code draws heavily from
 * that of R package 'bnstruct'
 * version 1.0.2
 * 
 * In Windows, it’s necessary that the Rtools 
 * executables directory 
 * (typically C:\Rtools\bin) and the C 
 * compiler executables directory (typically 
 * C:\Rtools\gcc-4.6.3\bin) are included in 
 * the Windows PATH environment variable. 
 * You may need to reboot Windows before R 
 * can recognise these values.
 */
#include <R.h>
#include <Rinternals.h>
#include <math.h>
// #include <errno.h>

/* =================================================================================== */
/* 
 * Input:
 * 'word': Parent set index in decimal.
 * 'bits':  Unsigned int dyn arr of length total no. of nodes in the data. 
 * 'size': Number of noeds in the data.
 * 
 * Output:
 * Returns 'count': No. of parents.
 * Modifies 'bits': Arr of node indices of parents. i^th elt of the arr
 * contains node index of the i^th parent. The first 'count' number of 
 * elts. contain node indices of parents.
 */
unsigned int get_bits( unsigned int word, unsigned int * bits, unsigned int size ) {
	unsigned int i, count = 0, bitmask = 1;
	for( i = 0; i < size; i++ ) {
	  
	  /* 
	   * Bit-wise binary AND, and
	   * bit-wise binary 'i' no. of left-shifts.
	   * Does not modify 'bitmask'.
	   */
	  if( word & bitmask<<i ) {
	    bits[count++] = i;
	  }
	}
	
	return count;
}

/* =================================================================================== */
/*
 * LL(G|D) - f(|D|)|B| 
 * |B| = \sum_{i=1}^{|V|} (r_i - 1) * q_i
 */
/*
 * Input parameters:
 * d = data.
 * n_nodes = number of nodes.
 * n_cases = number of replicates per node per time point.
 * ns = an array containing node sizes i.e. the discrete levels of nodes; 
 * 	ns[i] = number of discrete levls of the i-th node.
 * ni = node index of the target node.
 * pa = an array containing the indices of the candidate parent nodes.
 * n_pa = number of candidate parent nodes.
 * penalty = sparsity penalty to be deducted from the computed 
 * log likelihood for the given scoing function.
*/
double log_likelihood ( unsigned int * d, unsigned int n_nodes, unsigned int n_cases, unsigned int * ns, 
                        unsigned int ni, unsigned int * pa, unsigned int n_pa, double penalty ) {
  int i, j, index, elmt, stride;
  double acc, logl;
  
  // utility vectors
  int cum_prod_sizes[n_pa+2]; 
  cum_prod_sizes[0] = 1;
  cum_prod_sizes[1] = ns[ni];

  /*
   * Time complexity of the following for loop 
   * = n_pa
   * = number of candidate parents of the given target node
   * = O(M_f).
  */
  for( i = 0; i < n_pa; i++ )
    cum_prod_sizes[i+2] = cum_prod_sizes[i+1] * ns[pa[i]];
  
  int strides[n_pa + 1];
  strides[0] = ni * n_cases;

  /*
   * Time complexity of the following for loop 
   * = n_pa
   * = number of candidate parents of the given target node
   * = O(M_f).
  */	
  for( i = 0; i < n_pa; i++ )
    strides[i+1] = pa[i] * n_cases;
  
  // utility variables
  int prod_sizes = cum_prod_sizes[n_pa+1];
  int prod_sizes_pa = prod_sizes / ns[ni];
  // int n_na = 0;
  
  // vector holding counts (LARGE)
  double * counts = calloc( prod_sizes, sizeof(double) ); 
  
  // compute counts, skipping NAs
  /*
   * Number of iterations for the following for loop
   * = n_cases
   * = S.
  */
  for( i = 0; i < n_cases; i++ ) {
    index = 0;

    // sum using strides
    /*
     * Time complexity of the following for loop 
     * = (n_pa + 1)
     * = (number of candidate parents of the given target node + 1)
     * = (O(M_f) + 1).
    */
    for( j = 0; j < n_pa + 1; j++ ) {
      elmt = d[ i + strides[j] ];
      if( elmt == NA_INTEGER )
        //{
        //	n_na++;
        break;
      //}
      index += (elmt - 1) * cum_prod_sizes[j];
      // Rprintf("after index+\n");
    }
    // check if NA encountered
    if( j < n_pa + 1 )
      continue;
    
    counts[index] += 1;	
    // Rprintf("i = %d,\tj = %d\n", i, j);
  }
  // Rprintf("after loop\n");
  
  // correct for NAs
  //for( i = 0; i < prod_sizes; i++ )
  //	counts[i] += n_na * (counts[i] + prior) / ( n_cases - n_na + alpha );
  
  // compute log likelihood
  logl = 0.0;
  
  /*
   * Number of iterations for the following for loop
   * = prod_sizes_pa
   * = product of the numbers of discrete levels of 
   * 	the candidate parent nodes.
  */
  for( i = 0; i < prod_sizes_pa; i++ )
  {
    stride = i * ns[ni];
    acc    = counts[stride];

    /*
     * Time complexity of the following for loop 
     * = (ns[ni] - 1)
     * = (Number of discrete levels of the target node - 1).
    */    
    for( j = 1; j < ns[ni]; j++ )
    {
      acc  += counts[stride + j];
    }
    
    acc = log(acc + 1);
    
    /*
     * Time complexity of the following for loop 
     * = (ns[ni] - 1)
     * = (Number of discrete levels of the target node - 1).
    */
    for ( j = 1 ; j < ns[ni] ; j++ )
    {
      logl += counts[stride + j] * (log(counts[stride + j] + 1) - acc);
    }
    
  }
  
  logl -= penalty * prod_sizes_pa * (ns[ni] - 1);
  
  free(counts);
  // Rprintf("after free");
  
  return logl;
}

/* =================================================================================== */
double score_node_1( int* data, int ncols_data, int nrows_data, int* node_sizes, 
                     int ni, int* pars, int length_pars, int func) {
  double score;
  
  switch (func) {
    /* 
     * Only BIC scoring function is allowed till now
     */
    /*
     case 0 : score = bdeu_score( data, ncols_data, nrows_data, node_sizes,
     ni, pars, length_pars, ess );
     break;
     
     case 1 : // score = log_likelihood( data, ncols_data, nrows_data, node_sizes,
     // ni, pars, length_pars, 0.5*log(nrows_data) );
     // break;
     printf("\n");
     
     */
    
  case 2: // Rprintf("before ll\n");
    score = log_likelihood( data, ncols_data, nrows_data, node_sizes,
                                   ni, pars, length_pars, 1.0 );
    // Rprintf("after ll\n");
    break;
  // default: Rprintf("Only supports BIC scoring function till now.\n");
  // exit(EXIT_FAILURE);
  }
  
  return score;
}

/* =================================================================================== */
/*
 * aflml <- .Call("all_fam_log_marg_lik_v2", data, node.sizes, scoring.func)
 */
SEXP all_fam_log_marg_lik_v2(SEXP data, SEXP node_sizes, SEXP func, 
                             SEXP tgt_node_idx, SEXP src_node_idx, 
                             SEXP num_src_nodes) {
  
  // unsigned int i,j,n_pa,pos;
  
  // Begin: Get inputs
  unsigned int * d = INTEGER(data);
  // Rprintf("d\n");
  
  unsigned int n_nodes = ncols(data);
  
  // No. of obs.
  unsigned int n_ex = nrows(data);	 
  
  unsigned int * ns = INTEGER(node_sizes);
  // Rprintf("ns\n");
  // Impossible family mask
  // unsigned int * ifm = INTEGER(imp_fam_mask);
  // unsigned int pow_nodes = ncols(imp_fam_mask);
  
  // double alpha = *(REAL(iss));
  
  /*
  * 0 for BDeu
  * 1 for AIC
  * 2 for BIC
  */
  unsigned int scoring_func = *INTEGER(func);
  // Rprintf("scoring_func\n");
  
  int tgt = *INTEGER(tgt_node_idx);
  // Rprintf("tgt\n");
  
  int * pa = INTEGER(src_node_idx);
  // Rprintf("pa\n");
  
  int n_pa = *INTEGER(num_src_nodes);
  
  // unsigned int * pa = (unsigned int *) R_alloc( n_nodes, sizeof(unsigned int) );
  // End: Get inputs
  
  // allocate and initialize output
  SEXP result = PROTECT(allocVector(REALSXP, 1));
  
  // PROTECT( result = allocMatrix(REALSXP, n_nodes, pow_nodes) );
  // double * aflml = REAL(result);
  // for( i = 0; i < n_nodes*pow_nodes; i++ )
  // aflml[i] = R_NegInf;
  
  // Rprintf("before score\n");
  
  // compute log likelihood
  // REAL(result)[0] = -12.2932;
  REAL(result)[0] = score_node_1(d, n_nodes, n_ex, ns, tgt, pa, n_pa, scoring_func);
  
  // Rprintf("after score\n");
  // Rprintf("%f\n", REAL(result)[0]);
  
  // compute log likelihood
  // for( i = 0; i < n_nodes; i++ )
  // for( j = 0; j < pow_nodes; j++ ) {
  // pos = j*n_nodes + i;
  // if( ifm[pos] ) {
  /* 
   * 'n_pa' = No. of parents of node 'i'.
   * 'pa' = Arr of node indices of parents in the j^th parent set. 
   * k^th elt of the arr contains node index of the k^th parent. 
   * The first 'n_pa' number of 
   * elts. contain node indices of parents.
   */
  // Rprintf("get bits\n");
  // n_pa = get_bits( j, pa, n_nodes );
  
  /*
   * Calc score when the i-th node is parented by the j-th parent set.
   * The j^th parent set is represented by 'pa'.
   */			  
  // Rprintf("log lik, node %d, n parents %d\n",i,n_pa);
  // aflml[pos] = score_node_1(d, n_nodes, n_ex, ns, i, pa, n_pa, scoring_func);
  // bdeu_score( d, n_nodes, n_ex, ns, i, pa, n_pa, alpha );
  // Rprintf("end\n");
  // }
  // }
  
  /*
   * UNPROTECT() takes a single integer argument, n, and unprotects the 
   * last n objects that were protected. The number of protects and 
   * unprotects must match. If not, R will warn about a 
   * “stack imbalance in .Call”
   */
  UNPROTECT(1);
  return result;	
}
