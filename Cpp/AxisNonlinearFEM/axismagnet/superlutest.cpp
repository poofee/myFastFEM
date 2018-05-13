/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required
approvals from U.S. Dept. of Energy)

All rights reserved.

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*
 * -- SuperLU MT routine (version 3.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * and Xerox Palo Alto Research Center.
 * September 10, 2007
 *
 * Purpose:
 * ========
 *
 * This example illustrates how to use PDGSSVX to solve systems repeatedly
 * with the same sparsity pattern of matrix A.
 * In this case, the column permutation vector perm_c is computed once.
 * The following data structures will be reused in the subsequent call to
 * PDGSSVX: perm_c, etree, colcnt_h, part_super_h.
 *
 */
#include "superlutest.h"
#include "slu_mt_ddefs.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <QDebug>


/* Eat up the rest of the current line */
int_t dDumpLine(FILE *fp)
{
    int_t c;
    while ((c = fgetc(fp)) != '\n') ;
    return 0;
}

int_t dParseIntFormat(char *buf, int_t *num, int_t *size)
{
    char *tmp;

    tmp = buf;
    while (*tmp++ != '(') ;
    *num = atoi(tmp);
    while (*tmp != 'I' && *tmp != 'i') ++tmp;
    ++tmp;
    *size = atoi(tmp);
    return 0;
}

int_t dParseFloatFormat(char *buf, int_t *num, int_t *size)
{
    char *tmp, *period;

    tmp = buf;
    while (*tmp++ != '(') ;
    *num = atoi(tmp); /*sscanf(tmp, "%d", num);*/
    while (*tmp != 'E' && *tmp != 'e' && *tmp != 'D' && *tmp != 'd'
       && *tmp != 'F' && *tmp != 'f') {
        /* May find kP before nE/nD/nF, like (1P6F13.6). In this case the
           num picked up refers to P, which should be skipped. */
        if (*tmp=='p' || *tmp=='P') {
           ++tmp;
           *num = atoi(tmp); /*sscanf(tmp, "%d", num);*/
        } else {
           ++tmp;
        }
    }
    ++tmp;
    period = tmp;
    while (*period != '.' && *period != ')') ++period ;
    *period = '\0';
    *size = atoi(tmp); /*sscanf(tmp, "%2d", size);*/

    return 0;
}

int_t dReadVector(FILE *fp, int_t n, int_t *where, int_t perline, int_t persize)
{
    int_t i, j, item;
    char tmp, buf[100];

    i = 0;
    while (i < n) {
    fgets(buf, 100, fp);    /* read a line at a time */
    for (j=0; j<perline && i<n; j++) {
        tmp = buf[(j+1)*persize];     /* save the char at that place */
        buf[(j+1)*persize] = 0;       /* null terminate */
        item = atoi(&buf[j*persize]);
        buf[(j+1)*persize] = tmp;     /* recover the char at that place */
        where[i++] = item - 1;
    }
    }

    return 0;
}

int_t dReadValues(FILE *fp, int_t n, double *destination, int_t perline, int_t persize)
{
    int_t i, j, k, s;
    char tmp, buf[100];

    i = 0;
    while (i < n) {
    fgets(buf, 100, fp);    /* read a line at a time */
    for (j=0; j<perline && i<n; j++) {
        tmp = buf[(j+1)*persize];     /* save the char at that place */
        buf[(j+1)*persize] = 0;       /* null terminate */
        s = j*persize;
        for (k = 0; k < persize; ++k) /* No D_ format in C */
        if ( buf[s+k] == 'D' || buf[s+k] == 'd' ) buf[s+k] = 'E';
        destination[i++] = atof(&buf[s]);
        buf[(j+1)*persize] = tmp;     /* recover the char at that place */
    }
    }

    return 0;
}



void
dreadhb(FILE*fp,int_t *nrow, int_t *ncol, int_t *nonz,
    double **nzval, int_t **rowind, int_t **colptr)
{
/*
 * Purpose
 * =======
 *
 * Read a DOUBLE PRECISION matrix stored in Harwell-Boeing format
 * as described below.
 *
 * Line 1 (A72,A8)
 *  	Col. 1 - 72   Title (TITLE)
 *	Col. 73 - 80  Key (KEY)
 *
 * Line 2 (5I14)
 * 	Col. 1 - 14   Total number of lines excluding header (TOTCRD)
 * 	Col. 15 - 28  Number of lines for pointers (PTRCRD)
 * 	Col. 29 - 42  Number of lines for row (or variable) indices (INDCRD)
 * 	Col. 43 - 56  Number of lines for numerical values (VALCRD)
 *	Col. 57 - 70  Number of lines for right-hand sides (RHSCRD)
 *                    (including starting guesses and solution vectors
 *		       if present)
 *           	      (zero indicates no right-hand side data is present)
 *
 * Line 3 (A3, 11X, 4I14)
 *   	Col. 1 - 3    Matrix type (see below) (MXTYPE)
 * 	Col. 15 - 28  Number of rows (or variables) (NROW)
 * 	Col. 29 - 42  Number of columns (or elements) (NCOL)
 *	Col. 43 - 56  Number of row (or variable) indices (NNZERO)
 *	              (equal to number of entries for assembled matrices)
 * 	Col. 57 - 70  Number of elemental matrix entries (NELTVL)
 *	              (zero in the case of assembled matrices)
 * Line 4 (2A16, 2A20)
 * 	Col. 1 - 16   Format for pointers (PTRFMT)
 *	Col. 17 - 32  Format for row (or variable) indices (INDFMT)
 *	Col. 33 - 52  Format for numerical values of coefficient matrix (VALFMT)
 * 	Col. 53 - 72 Format for numerical values of right-hand sides (RHSFMT)
 *
 * Line 5 (A3, 11X, 2I14) Only present if there are right-hand sides present
 *    	Col. 1 	      Right-hand side type:
 *	         	  F for full storage or M for same format as matrix
 *    	Col. 2        G if a starting vector(s) (Guess) is supplied. (RHSTYP)
 *    	Col. 3        X if an exact solution vector(s) is supplied.
 *	Col. 15 - 28  Number of right-hand sides (NRHS)
 *	Col. 29 - 42  Number of row indices (NRHSIX)
 *          	      (ignored in case of unassembled matrices)
 *
 * The three character type field on line 3 describes the matrix type.
 * The following table lists the permitted values for each of the three
 * characters. As an example of the type field, RSA denotes that the matrix
 * is real, symmetric, and assembled.
 *
 * First Character:
 *	R Real matrix
 *	C Complex matrix
 *	P Pattern only (no numerical values supplied)
 *
 * Second Character:
 *	S Symmetric
 *	U Unsymmetric
 *	H Hermitian
 *	Z Skew symmetric
 *	R Rectangular
 *
 * Third Character:
 *	A Assembled
 *	E Elemental matrices (unassembled)
 *
 */

    int_t i, numer_lines, rhscrd = 0;
    int_t tmp, colnum, colsize, rownum, rowsize, valnum, valsize;
    char buf[100], type[4], key[10];
    //FILE *fp;

    //fp = stdin;

    /* Line 1 */
    fscanf(fp, "%72c", buf); buf[72] = 0;
    printf("Title: %s", buf);
    fscanf(fp, "%8c", key);  key[8] = 0;
    printf("Key: %s\n", key);
    dDumpLine(fp);

    /* Line 2 */
    for (i=0; i<5; i++) {
    fscanf(fp, "%14c", buf); buf[14] = 0;
    tmp = atoi(buf); /*sscanf(buf, "%d", &tmp);*/
    if (i == 3) numer_lines = tmp;
    if (i == 4 && tmp) rhscrd = tmp;
    }
    dDumpLine(fp);

    /* Line 3 */
    fscanf(fp, "%3c", type);
    fscanf(fp, "%11c", buf); /* pad */
    type[3] = 0;
#if ( DEBUGlevel>=1 )
    printf("Matrix type %s\n", type);
#endif

    fscanf(fp, "%14c", buf); *nrow = atoi(buf);
    fscanf(fp, "%14c", buf); *ncol = atoi(buf);
    fscanf(fp, "%14c", buf); *nonz = atoi(buf);
    fscanf(fp, "%14c", buf); tmp = atoi(buf);

    if (tmp != 0)
      printf("This is not an assembled matrix!\n");
    if (*nrow != *ncol)
    printf("Matrix is not square.\n");
    dDumpLine(fp);

    /* Line 4: format statement */
    fscanf(fp, "%16c", buf);
    dParseIntFormat(buf, &colnum, &colsize);
    fscanf(fp, "%16c", buf);
    dParseIntFormat(buf, &rownum, &rowsize);
    fscanf(fp, "%20c", buf);
    dParseFloatFormat(buf, &valnum, &valsize);
    fscanf(fp, "%20c", buf);
    dDumpLine(fp);

    /* Line 5: right-hand side */
    if ( rhscrd ) dDumpLine(fp); /* skip RHSFMT */

#if ( DEBUGlevel>=1 )
    printf("%d rows, %d nonzeros\n", *nrow, *nonz);
    printf("colnum %d, colsize %d\n", colnum, colsize);
    printf("rownum %d, rowsize %d\n", rownum, rowsize);
    printf("valnum %d, valsize %d\n", valnum, valsize);
#endif

    /* Allocate storage for the three arrays ( nzval, rowind, colptr ) */
    dallocateA(*ncol, *nonz, nzval, rowind, colptr);

    dReadVector(fp, *ncol+1, *colptr, colnum, colsize);
    dReadVector(fp, *nonz, *rowind, rownum, rowsize);
    if ( numer_lines ) {
        dReadValues(fp, *nonz, *nzval, valnum, valsize);
    }


}
void superlumttest()
{
    //修改为自己目录下的文件
    FILE * fp = fopen("/Users/poofee/big.rua","r");
    if(!fp){
        printf("Error open file!");
        return;
    }else{
        qDebug()<<"Opened file successfully!";
    }

    SuperMatrix A, A1, L, U;
    SuperMatrix B, B1, X;
    NCformat    *Astore;
    SCPformat   *Lstore;
    NCPformat   *Ustore;
    int_t         nprocs;
    fact_t      fact;
    trans_t     trans;
    yes_no_t    refact, usepr;
    equed_t     equed;
    double      *a, *a1;
    int_t         *asub, *xa, *asub1, *xa1;
    int_t         *perm_c; /* column permutation vector */
    int_t         *perm_r; /* row permutations from partial pivoting */
    void        *work;
    superlumt_options_t superlumt_options;
    int_t         info, lwork, nrhs, ldx, panel_size, relax;
    int_t         m, n, nnz, permc_spec, i;
    double      *rhsb, *rhsx, *xact;
    double      *R, *C;
    double      *ferr, *berr;
    double      u, drop_tol, rpg, rcond;
    superlu_memusage_t superlu_memusage;
    //void parse_command_line();

    /* Default parameters to control factorization. */
    nprocs = 1;
    fact  = EQUILIBRATE;
    trans = NOTRANS;
    equed = NOEQUIL;
    refact= NO;
    panel_size = sp_ienv(1);
    relax = sp_ienv(2);
    u     = 1.0;
    usepr = NO;
    drop_tol = 0.0;
    lwork = 0;
    nrhs  = 1;

    /* Command line options to modify default behavior. */
    ///parse_command_line(argc, argv, &nprocs, &lwork, &panel_size, &relax,
    //	       &u, &fact, &trans, &refact, &equed);

    if ( lwork > 0 ) {
    work = SUPERLU_MALLOC(lwork);
    printf("Use work space of size LWORK = " IFMT " bytes\n", lwork);
    if ( !work ) {
        SUPERLU_ABORT("DLINSOLX: cannot allocate work[]");
    }
    }

#if ( PRNTlevel==1 )
    cpp_defs();
    printf("int_t %d bytes\n", sizeof(int_t));
#endif

#define HB
#if defined( DEN )
    m = n;
    nnz = n * n;
    dband(n, n, nnz, &a, &asub, &xa);
#elif defined( BAND )
    m = n;
    nnz = (2*b+1) * n;
    dband(n, b, nnz, &a, &asub, &xa);
#elif defined( BD )
    nb = 5;
    bs = 200;
    m = n = bs * nb;
    nnz = bs * bs * nb;
    dblockdiag(nb, bs, nnz, &a, &asub, &xa);
#elif defined( HB )
    dreadhb(fp,&m, &n, &nnz, &a, &asub, &xa);
#else
    dreadtriple(&m, &n, &nnz, &a, &asub, &xa);
#endif

    fclose(fp);
    qDebug()<<sizeof(int)<<sizeof(int_t);

    if ( !(a1 = doubleMalloc(nnz)) ) qDebug()<<"Malloc fails for a1[].";
    if ( !(asub1 = intMalloc(nnz)) ) qDebug()<<"Malloc fails for asub1[].";
    if ( !(xa1 = intMalloc(n+1)) ) qDebug()<<"Malloc fails for xa1[].";
    for (i = 0; i < nnz; ++i) {
        a1[i] = a[i];
    asub1[i] = asub[i];
    }
    for (i = 0; i < n+1; ++i) xa1[i] = xa[i];

    dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
    Astore = (NCformat*)A.Store;
    qDebug()<<Astore->nnz;
    printf("Dimension " IFMT "x" IFMT "; # nonzeros " IFMT "\n", A.nrow, A.ncol, Astore->nnz);

    if (!(rhsb = doubleMalloc(m * nrhs))) SUPERLU_ABORT("Malloc fails for rhsb[].");
    if (!(rhsx = doubleMalloc(m * nrhs))) SUPERLU_ABORT("Malloc fails for rhsx[].");
    dCreate_Dense_Matrix(&B, m, nrhs, rhsb, m, SLU_DN, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&X, m, nrhs, rhsx, m, SLU_DN, SLU_D, SLU_GE);
    xact = doubleMalloc(n * nrhs);
    ldx = n;
    dGenXtrue(n, nrhs, xact, ldx);
    dFillRHS(trans, nrhs, xact, ldx, &A, &B);

    if (!(perm_r = intMalloc(m))) SUPERLU_ABORT("Malloc fails for perm_r[].");
    if (!(perm_c = intMalloc(n))) SUPERLU_ABORT("Malloc fails for perm_c[].");
    if (!(R = (double *) SUPERLU_MALLOC(A.nrow * sizeof(double))))
        SUPERLU_ABORT("SUPERLU_MALLOC fails for R[].");
    if ( !(C = (double *) SUPERLU_MALLOC(A.ncol * sizeof(double))) )
        SUPERLU_ABORT("SUPERLU_MALLOC fails for C[].");
    if ( !(ferr = (double *) SUPERLU_MALLOC(nrhs * sizeof(double))) )
        SUPERLU_ABORT("SUPERLU_MALLOC fails for ferr[].");
    if ( !(berr = (double *) SUPERLU_MALLOC(nrhs * sizeof(double))) )
        SUPERLU_ABORT("SUPERLU_MALLOC fails for berr[].");

    /*
     * Get column permutation vector perm_c[], according to permc_spec:
     *   permc_spec = 0: natural ordering
     *   permc_spec = 1: minimum degree ordering on structure of A'*A
     *   permc_spec = 2: minimum degree ordering on structure of A'+A
     *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
     */
    permc_spec = 1;
    get_perm_c(permc_spec, &A, perm_c);

    superlumt_options.nprocs = nprocs;
    superlumt_options.fact = fact;
    superlumt_options.trans = trans;
    superlumt_options.refact = refact;
    superlumt_options.panel_size = panel_size;
    superlumt_options.relax = relax;
    superlumt_options.usepr = usepr;
    superlumt_options.drop_tol = drop_tol;
    superlumt_options.diag_pivot_thresh = u;
    superlumt_options.SymmetricMode = NO;
    superlumt_options.PrintStat = NO;
    superlumt_options.perm_c = perm_c;
    superlumt_options.perm_r = perm_r;
    superlumt_options.work = work;
    superlumt_options.lwork = lwork;
    if ( !(superlumt_options.etree = intMalloc(n)) )
    SUPERLU_ABORT("Malloc fails for etree[].");
    if ( !(superlumt_options.colcnt_h = intMalloc(n)) )
    SUPERLU_ABORT("Malloc fails for colcnt_h[].");
    if ( !(superlumt_options.part_super_h = intMalloc(n)) )
    SUPERLU_ABORT("Malloc fails for colcnt_h[].");

    /* ------------------------------------------------------------
       WE SOLVE THE LINEAR SYSTEM FOR THE FIRST TIME: AX = B
       ------------------------------------------------------------*/
    pdgssvx(nprocs, &superlumt_options, &A, perm_c, perm_r,
        &equed, R, C, &L, &U, &B, &X, &rpg, &rcond,
        ferr, berr, &superlu_memusage, &info);

    if ( info == 0 || info == n+1 ) {

    printf("Recip. pivot growth = %e\n", rpg);
    printf("Recip. condition number = %e\n", rcond);
    printf("%8s%16s%16s\n", "rhs", "FERR", "BERR");
    for (i = 0; i < nrhs; ++i) {
        printf(IFMT "%16e%16e\n", i+1, ferr[i], berr[i]);
    }

        Lstore = (SCPformat *) L.Store;
        Ustore = (NCPformat *) U.Store;
    printf("No of nonzeros in factor L = " IFMT "\n", Lstore->nnz);
        printf("No of nonzeros in factor U = " IFMT "\n", Ustore->nnz);
        printf("No of nonzeros in L+U = " IFMT "\n", Lstore->nnz + Ustore->nnz - n);
    printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions " IFMT "\n",
           superlu_memusage.for_lu/1e6, superlu_memusage.total_needed/1e6,
           superlu_memusage.expansions);

    fflush(stdout);

    } else if ( info > 0 && lwork == -1 ) {
        printf("** Estimated memory: " IFMT " bytes\n", info - n);
    }

    printf("First system: pdgssvx(): info " IFMT "\n----\n", info);

    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);

    /* ------------------------------------------------------------
       NOW WE SOLVE ANOTHER LINEAR SYSTEM: A1*X = B1
       ONLY THE SPARSITY PATTERN OF A1 IS THE SAME AS THAT OF A.
       ------------------------------------------------------------*/
    superlumt_options.refact = YES;
    dCreate_CompCol_Matrix(&A1, m, n, nnz, a1, asub1, xa1, SLU_NC, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&B1, m, nrhs, rhsb, m, SLU_DN, SLU_D, SLU_GE);

    pdgssvx(nprocs, &superlumt_options, &A1, perm_c, perm_r,
        &equed, R, C, &L, &U, &B1, &X, &rpg, &rcond,
        ferr, berr, &superlu_memusage, &info);

    if ( info == 0 || info == n+1 ) {

    printf("Recip. pivot growth = %e\n", rpg);
    printf("Recip. condition number = %e\n", rcond);
    printf("%8s%16s%16s\n", "rhs", "FERR", "BERR");
    for (i = 0; i < nrhs; ++i) {
        printf(IFMT "%16e%16e\n", i+1, ferr[i], berr[i]);
    }

        Lstore = (SCPformat *) L.Store;
        Ustore = (NCPformat *) U.Store;
    printf("No of nonzeros in factor L = " IFMT "\n", Lstore->nnz);
        printf("No of nonzeros in factor U = " IFMT "\n", Ustore->nnz);
        printf("No of nonzeros in L+U = " IFMT "\n", Lstore->nnz + Ustore->nnz - n);
    printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions " IFMT "\n",
           superlu_memusage.for_lu/1e6, superlu_memusage.total_needed/1e6,
           superlu_memusage.expansions);

    fflush(stdout);

    } else if ( info > 0 && lwork == -1 ) {
        printf("** Estimated memory: " IFMT " bytes\n", info - n);
    }

    printf("Second system: pdgssvx(): info " IFMT "\n", info);

    SUPERLU_FREE (rhsb);
    SUPERLU_FREE (rhsx);
    SUPERLU_FREE (xact);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    SUPERLU_FREE (R);
    SUPERLU_FREE (C);
    SUPERLU_FREE (ferr);
    SUPERLU_FREE (berr);
    Destroy_CompCol_Matrix(&A1);
    Destroy_SuperMatrix_Store(&B1);
    Destroy_SuperMatrix_Store(&X);
    SUPERLU_FREE (superlumt_options.etree);
    SUPERLU_FREE (superlumt_options.colcnt_h);
    SUPERLU_FREE (superlumt_options.part_super_h);
    if ( lwork == 0 ) {
        Destroy_SuperNode_SCP(&L);
        Destroy_CompCol_NCP(&U);
    } else if ( lwork > 0 ) {
        SUPERLU_FREE(work);
    }
}


