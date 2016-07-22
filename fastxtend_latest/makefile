all : fastx_duplicatedReads fastx_clean fastx_mergepairs fastx_stats

fastx_duplicatedReads : libfastx_comp
	cd fastx_duplicatedReads; $(MAKE) ;
		
fastx_clean : libfastx_comp
	cd fastx_clean; $(MAKE) ;	
	
fastx_mergepairs : libfastx_comp
	cd fastx_mergepairs; $(MAKE) ;
	
fastx_stats : libfastx_comp
	cd fastx_stats; $(MAKE) ;
		
libfastx_comp :
	cd libfastx; $(MAKE) ;
	
install :
	cd fastx_duplicatedReads; $(MAKE) install ;
	cd fastx_clean; $(MAKE) install ;
	cd fastx_mergepairs; $(MAKE) install ;
	cd fastx_stats; $(MAKE) install ;

clean :
	cd fastx_duplicatedReads; $(MAKE) clean ;
	cd fastx_clean; $(MAKE) clean ;
	cd fastx_mergepairs; $(MAKE) clean ;
	cd fastx_stats; $(MAKE) clean ;
	cd libfastx; $(MAKE) clean ;
