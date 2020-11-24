bin_PROGRAMS=pulsarfold

AM_CPPFLAGS=-I$(top_srcdir)/include
LDFLAGS=-L$(top_srcdir)/src/container -L$(top_srcdir)/src/formats -L$(top_srcdir)/src/utils -L$(top_srcdir)/src/module
LDADD=-lboost_program_options -lformats -lutils -lcontainer -lmodule

pulsarfold_SOURCES=pulsarfold.cpp pulsarfolder.cpp predictor.cpp
