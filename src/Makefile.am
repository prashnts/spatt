# $Id: Makefile.am 553 2005-11-24 13:11:19Z gnuel $

bin_PROGRAMS = sspatt ldspatt oldxspatt cpspatt gspatt xspatt aspatt

noinst_LTLIBRARIES = libspatt.la

libspatt_la_SOURCES = \
                 spattglobals.h \
		 spattparameters.h spattparameters.cc \
                 spattparameterparser.h spattparameterparser.cc \
		 alphabet.h alphabet.cc\
		 sequence.h sequence.cc\
		 count.h count.cc\
		 word.h word.cc\
		 input.h input.cc\
		 pattern.h pattern.cc\
		 markov.h markov.cc\
		 cdf.h cdf.c \
                 cp.h cp.cc \
		 fortran.h fortran.cc\
		 process.h process.cc\
		 stat.h stat.cc \
		 sstat.h sstat.cc

sspatt_SOURCES = \
		sspattparameters.h sspattparameters.cc \
		sspattparameterparser.h sspattparameterparser.cc \
		sspatt.cc

aspatt_SOURCES = \
		amarkov.h amarkov.cc \
		fortran.h fortran.cc \
		cdf.h cdf.c \
		timer.h \
		cp.h cp.cc \
		aspatt.cc

oldxspatt_SOURCES = \
		xspattparameters.h xspattparameters.cc \
		xspattparameterparser.h xspattparameterparser.cc \
		cpstat.h cpstat.cc \
		xstat.h xstat.cc \
		x_WordFam.h x_WordFam.cc \
		x_PSucceed.h x_PSucceed.cc \
		x_PAppearFast.h x_PAppearFast.cc \
		oldxspatt.cc

xspatt_SOURCES = \
		xspattparameters.h xspattparameters.cc \
		xspattparameterparser.h xspattparameterparser.cc \
		transition.h transition.cc \
		newxstat.h newxstat.cc \
	       	xspatt.cc

ldspatt_SOURCES = \
		ldspattparameters.h ldspattparameters.cc \
		ldspattparameterparser.h ldspattparameterparser.cc \
		transition.h transition.cc \
		ldstat.h ldstat.cc \
	       	ldspatt.cc

cpspatt_SOURCES = \
                 cpspattparameters.h cpspattparameters.cc \
                 cpspattparameterparser.h cpspattparameterparser.cc \
		 cpstat.h cpstat.cc \
                 cpspatt.cc

gspatt_SOURCES = \
                 gspattparameters.h gspattparameters.cc \
                 gspattparameterparser.h gspattparameterparser.cc \
		 gstat.h gstat.cc \
                 gspatt.cc


LDADD = libspatt.la -L../arpack -lArpack -L../lapack -lLapack -L../blas -lBlas @FLIBS@ @GSL_LIBS@

INCLUDES = -I../tclap-1.0.5/include @GSL_CFLAGS@
