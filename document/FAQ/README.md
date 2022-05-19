## This is a SAW FAQ | Frequently asked questions about running SAW

## FAQ
**Question:** When will cellBin be available on SAW?

**Answer:** 

##
**Q:** no spaces in bind mount specification

**A:** 
In the SAW pipeline scripts for single lane analysis at [1], there is a bug. The SINGULARITY_BIND environment variable (a comma-separated list) is set with spaces in front of two arguments. This causes the bind mounting to fail.
This bug is only present in the single-lane analysis script - the multi-lane analysis script does not have this issue.

##
**Q:** 

**A:**
##
