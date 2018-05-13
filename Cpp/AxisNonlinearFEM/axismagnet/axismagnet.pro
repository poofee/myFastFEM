#-------------------------------------------------
#
# Project created by QtCreator 2018-03-16T15:27:20
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = axismagnet
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    axismagnet2d.cpp \
    datatype.cpp \
    superlutest.cpp

HEADERS  += mainwindow.h \
    axismagnet2d.h \
    datatype.h \
    superlutest.h

FORMS    += mainwindow.ui

macx {
    message("The OS is macOS.")
#librarys
LIBS += \
    -L../axismagnet/SuperLU_MT_3.1/lib -lsuperlu_mt_PTHREAD \
    -L../axismagnet/OpenBLAS/lib/lib -lopenblas \
    -lpthread \

#additional include path
INCLUDEPATH += \
    ./SuperLU_MT_3.1/SRC \
    ./OpenBLAS/include
}
