# CTPOL_MD

MD script to run CTPOL model on OpenMM.

## Required files:

* Original OPLSAA parameter file, must be named `OPLS-AA.xml`.
* Input pdb, either as 2nd [optional] argument to script, or a file called `input.pdb`.

### For CTPOL,
* Modified parameter file, 1st argument to script, must have `ffaffurr-oplsaa` in name.
* A parameter file for CT and Pol, must be named `CustomForce.xml`

# Usage:

`python ffaffur-md-1step.py <parameter-file> [pdb]`

