begin program.
import spss, spssaux
reload(ETTEST)
end program.

ETTEST GROUPS=iv 
   /GROUP GONE=2 GTWO=1
   /WITH DVLIST=dv otherdv
   /OPTIONS MISSING=analysis CI=90.
