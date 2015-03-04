
!! careful, ECSFLAG should only be defined if REALGO isn't

!! REALGO determines DATATYPE - type of most of data in calculation.

!! CNORMFLAG determines CONJUGATE, determines dot product.
!! goes into mean fields and denmat, and dot, and therefore gramschmidt
!! goes with EIGEN (for denmat)  CONJUGATE/REALGO define CNORMTYPE
!!  which is type of denvals (only).
!! 

!! ECSFLAG defines DATAECS.
!!  mostly used deep in program (PSC / MYGETRF)
!!  but determines CONFIGEIG; DATAECS is type of tempconfigvals
!!  and determines ECSEIG for spfs
!!  except for myprojectmod/TWOE.f90/PSC.f90, DATAECS used only sparingly
!!  not used in any routines with MYGEM* - OK!


#define REALxxGO
#define ECSFLAG
#define CNOxxRMFLAG

#include "Definitions.ALL"



