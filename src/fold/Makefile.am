bin_PROGRAMS=pulsarfold pulsarfold_fil

AM_CPPFLAGS=-I$(top_srcdir)/include
LDFLAGS=-L$(top_srcdir)/src/container -L$(top_srcdir)/src/formats -L$(top_srcdir)/src/utils -L$(top_srcdir)/src/module
LDADD=-lboost_program_options -lmodule -lformats -lutils -lcontainer

pulsarfold_SOURCES=pulsarfold.cpp pulsarfolder.cpp predictor.cpp
pulsarfold_fil_SOURCES=pulsarfold_fil.cpp pulsarfolder.cpp predictor.cpp
