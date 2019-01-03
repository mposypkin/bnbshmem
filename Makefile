all dep clean indent tests::
	cd paralsolver && $(MAKE) $@ && cd .. 
	cd bnbatomic && $(MAKE) $@ && cd ..
	cd bnbomp && $(MAKE) $@ && cd ..
	cd lpamigo && $(MAKE) $@ && cd ..

doc: indent doxy

doxy:
	mkdir -p doc/html &&\
	doxygen doxy.conf

clean::
	rm -rf *~ PI* core bin/* obj/* tmp *.log