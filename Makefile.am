# $Id: Makefile.am 438 2005-07-20 12:09:50Z mhoebeke $
AUTOMAKE_OPTIONS = gnu
SUBDIRS = doc arpack blas lapack src tests
EXTRA_DIST = spatt.spec tclap-1.0.5
CLEANFILES = *~

STATICLDFLAGS=-static $(LDFLAGS)

STATICPROGRAMS=`find src -type f -perm +0111`

STATICFILES=\
        AUTHORS INSTALL README TODO BUGS\
        doc/man/*.1 \
	$(STATICPROGRAMS)

STATICDIR=@PACKAGE@-@VERSION@_static_$(host_triplet)


dist-hook :
	find $(distdir) -type d -name .svn | xargs rm -rf

static :
	rm -rf $(STATICDIR)
	rm -rf $(STATICPROGRAMS)
	make LDFLAGS=$(STATICLDFLAGS)
	mkdir $(STATICDIR)
	strip $(STATICPROGRAMS)
	@for file in $(STATICFILES) ;\
		do cp $$file $(STATICDIR) ;\
	done ;
	tar zcvf $(STATICDIR).tgz $(STATICDIR)
	rm -rf $(STATICDIR)
	rm -rf $(STATICPROGRAMS)
