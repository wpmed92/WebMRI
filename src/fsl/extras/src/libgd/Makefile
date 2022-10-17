include ${FSLCONFDIR}/default.mk

PROJNAME = libgd

OPTFLAGS = -O3 -DHAVE_LIBPNG

USRINCFLAGS = -I${INC_PNG} -I${INC_ZLIB}

TESTXFILES=gdtest gddemo gd2time gdtestft testac pngtogd pngtogd2 gdtopng gd2topng gd2copypal gdparttopng

OBJS = gd.o gd_gd.o gd_gd2.o gd_io.o gd_io_dp.o gd_io_file.o gd_ss.o gd_io_ss.o gd_png.o gdxpm.o gdfontt.o gdfonts.o gdfontmb.o gdfontl.o gdfontg.o gdtables.o gdft.o gdkanji.o wbmp.o gd_wbmp.o gdhelpers.o gd_topal.o

LIBS = -lpng -lz -lm

all: libgd.a 

test: $(TESTXFILES)

libgd.a: ${OBJS}
	${AR} -r libgd.a ${OBJS}


gddemo: gddemo.o
	$(CC)  ${CXXFLAGS} ${LDFLAGS} -o $@ gddemo.o ${LIBS} -lgd

testac: testac.o
	$(CC)  ${CXXFLAGS} ${LDFLAGS} -o $@ testac.o ${LIBS} -lgd

pngtogd: pngtogd.o
	$(CC)  ${CXXFLAGS} ${LDFLAGS} -o $@ pngtogd.o ${LIBS} -lgd

webpng: webpng.o
	$(CC)  ${CXXFLAGS} ${LDFLAGS} -o $@ webpng.o ${LIBS} -lgd
 
pngtogd2: pngtogd2.o
	$(CC)  ${CXXFLAGS} ${LDFLAGS} -o $@ pngtogd2.o ${LIBS} -lgd

gdtopng: gdtopng.o
	$(CC)  ${CXXFLAGS} ${LDFLAGS} -o $@ gdtopng.o ${LIBS} -lgd
 
gd2topng: gd2topng.o
	$(CC)  ${CXXFLAGS} ${LDFLAGS} -o $@ gd2topng.o ${LIBS} -lgd
 
gd2copypal: gd2copypal.o
	$(CC)  ${CXXFLAGS} ${LDFLAGS} -o $@ gd2copypal.o ${LIBS} -lgd

gdparttopng: gdparttopng.o
	$(CC)  ${CXXFLAGS} ${LDFLAGS} -o $@ gdparttopng.o ${LIBS} -lgd
 
gdtest: gdtest.o
	$(CC)  ${CXXFLAGS} ${LDFLAGS} -o $@ gdtest.o ${LIBS} -lgd
 
gd2time: gd2time.o
	$(CC)  ${CXXFLAGS} ${LDFLAGS} -o $@ gd2time.o ${LIBS} -lgd

gdtestft: gdtestft.o
	$(CC)  ${CXXFLAGS} ${LDFLAGS} -o $@ gdtestft.o ${LIBS} -lgd

