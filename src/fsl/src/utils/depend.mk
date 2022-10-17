check.o: check.cc options.h
functions.o: functions.cc options.h
log.o: log.cc log.h \
  /Users/ahmedharmouche/Documents/fsl/extras/include/newmat/newmatap.h \
  /Users/ahmedharmouche/Documents/fsl/extras/include/newmat/newmat.h \
  /Users/ahmedharmouche/Documents/fsl/extras/include/newmat/include.h \
  /Users/ahmedharmouche/Documents/fsl/extras/include/newmat/boolean.h \
  /Users/ahmedharmouche/Documents/fsl/extras/include/newmat/myexcept.h \
  /Users/ahmedharmouche/Documents/fsl/extras/include/newmat/newmatio.h
matches.o: matches.cc options.h
options.o: options.cc options.h
opttst.o: opttst.cc options.h
parse.o: parse.cc options.h
time_tracer.o: time_tracer.cc time_tracer.h
usage.o: usage.cc options.h
