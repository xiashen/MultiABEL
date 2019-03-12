#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/*
 The following name(s) appear with different usages
 e.g., with different numbers of arguments:
 
 MultiSummaryLoopPrecise
 
 This needs to be resolved in the tables and any declarations.
*/

/* FIXME:
 Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(multisummaryloop)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
				      extern void F77_NAME(multisummaryloopdirect)(void *, void *, void *, void *, void *, void *, void *);
										  extern void F77_NAME(multisummaryloopinbred)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
															      extern void F77_NAME(multisummaryloopprecise)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

																					   static const R_FortranMethodDef FortranEntries[] = {
																					     {"multisummaryloop",        (DL_FUNC) &F77_NAME(multisummaryloop),        16},
																					     {"multisummaryloopdirect",  (DL_FUNC) &F77_NAME(multisummaryloopdirect),   7},
																					     {"multisummaryloopinbred",  (DL_FUNC) &F77_NAME(multisummaryloopinbred),  16},
																					     {"multisummaryloopprecise", (DL_FUNC) &F77_NAME(multisummaryloopprecise), 16},
																					     {NULL, NULL, 0}
																					   };

void R_init_MultiABEL(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
