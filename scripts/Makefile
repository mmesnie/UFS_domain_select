URL=https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/

BASEDIR=~/ufs_DOMAIN_select/build/DATA-2.2.0/input_model_data/FV3GFS/grib2/${DATE}${CYCLE}

all: check download

check:
	@if [ "${DATE}" = "" ] || [ "${CYCLE}" = "" ] || [ ${LEN} = "" ] || [ "${LBC}" = "" ] ; then \
		echo "usage: make DATE=<yyyymmdd> CYCLE=<00|06|12|18> LEN=<int> LBC=<int>" ; \
		exit 1; \
	fi

basedir:
	@if [ ! -d $(BASEDIR) ]; then \
		mkdir -p $(BASEDIR); \
		cd $(BASEDIR)/..; \
		echo created $(BASEDIR); \
	else \
		echo already created $(BASEDIR); \
	fi

download: basedir
	@cd $(BASEDIR); \
	files=`seq 0 $(LBC) $(LEN) | \
	      awk '{if ($$1<10) printf("00%i\n", $$1); else printf("0%i\n", $$1)}'`; \
	for i in $$files; do \
		file=gfs.t$(CYCLE)z.pgrb2.0p25.f$$i; \
		if [ ! -f $$file ]; then \
			wget $(URL)/gfs.$(DATE)/$(CYCLE)/atmos/$$file; \
		else \
			echo already downloaded $$file; \
		fi; \
	done
