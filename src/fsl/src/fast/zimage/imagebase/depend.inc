image.o: image.cpp image.h imagebase.h ../include/mydefine.h pixfmt.h \
 rgbpixels.h ../include/common.h ../include/../common/commfunc.h \
 ../include/../common/point.h bitmap.h palette.h iterator.h \
 ../math/mathcommon.h ../math/matrix.h ../math/vector.h \
 ../math/vectorD.h ../math/iterator.h ../math/matrix2D.h \
 ../math/matrix3D.h ../math/matrix4D.h imagehelper.h
imagebase.o: imagebase.cpp imagebase.h ../include/mydefine.h pixfmt.h \
 rgbpixels.h ../include/common.h ../include/../common/commfunc.h \
 ../include/../common/point.h bitmap.h palette.h
imagehelper.o: imagehelper.cpp imagehelper.h image.h imagebase.h \
 ../include/mydefine.h pixfmt.h rgbpixels.h ../include/common.h \
 ../include/../common/commfunc.h ../include/../common/point.h bitmap.h \
 palette.h iterator.h ../math/mathcommon.h ../math/matrix.h \
 ../math/vector.h ../math/vectorD.h ../math/iterator.h \
 ../math/matrix2D.h ../math/matrix3D.h ../math/matrix4D.h
palette.o: palette.cpp palette.h ../include/mydefine.h bitmap.h
