TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    library.cpp
LIBS += -L/usr/local/lib
LIBS += -larmadillo -llapack -lblas
INCLUDEPATH += /usr/local/include

HEADERS += \
    library.h
