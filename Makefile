ECM_OPTS=
MSIEVE_OPTS=
GGNFS_OPTS=
YAFU_OPTS=
GMP_OPTS=

MAKE_OPTS=

PKG_DIR=$(shell cd '$( dirname '${BASH_SOURCE[0]}' )' && pwd)
PREFIX=$(PKG_DIR)/prefix
export PKG_DIR PREFIX

GMP_DIR=$(PKG_DIR)/gmp-6.0.0
ECM_DIR=$(PKG_DIR)/ecm-6.4.4
MSIEVE_DIR=$(PKG_DIR)/msieve-1.52
GGNFS_DIR=$(PKG_DIR)/ggnfs
YAFU_DIR=$(PKG_DIR)/yafu-1.34.3

LASIEVE_DIR=$(PKG_DIR)/lasieve_bin

export GMP_DIR ECM_DIR MSIEVE_DIR GGNFS_DIR YAFU_DIR LASIEVE_DIR

PKG_CFLAGS=-L$(PREFIX)/lib -I$(PREFIX)/include
PKG_LDFLAGS=-L$(PREFIX)/lib
export PKG_CFLAGS PKG_LDFLAGS
CFLAGS=$(PKG_CFLAGS)
LDFLAGS=$(PKG_LDFLAGS)
#export CFLAGS LDFLAGS

# Cleaning Tasks

clean_gmp:
	-$(MAKE) -C $(GMP_DIR) uninstall clean

clean_ecm:
	-$(MAKE) -C $(ECM_DIR) uninstall clean

clean_msieve:
	$(MAKE) -C $(MSIEVE_DIR) clean

clean_ggnfs:
	$(MAKE) -C $(GGNFS_DIR) clean

clean_yafu:
	$(MAKE) -C $(YAFU_DIR) clean
	rm -f $(PREFIX)/bin/yafu

clean: clean_gmp clean_ecm clean_msieve clean_ggnfs clean_yafu

# Subtasks

populate_prefix:
	
	mkdir -p $(PREFIX)/include $(PREFIX)/lib $(PREFIX)/bin $(PREFIX)/share

gmp: populate_prefix
	cd $(GMP_DIR) && \
	$(GMP_DIR)/configure --prefix=$(PREFIX) $(GMP_OPTS) && \
	$(MAKE) $(MAKE_OPTS) -C $(GMP_DIR) install

ecm: gmp
	cd $(ECM_DIR) && \
	$(ECM_DIR)/configure --prefix=$(PREFIX) $(EMC_OPTS) && \
	$(MAKE) $(MAKE_OPTS) -C $(ECM_DIR) install

msieve: ecm
	cd $(MSIEVE_DIR) && \
	$(MAKE) $(MAKE_OPTS) -C $(MSIEVE_DIR) all NO_ZLIB=1 ECM=1 $(MSIEVE_OPTS)

lasieve_binaries:
	-mkdir $(GGNFS_DIR)/bin
	cp $(LASIEVE_DIR)/gnfs-lasieve4I1* $(GGNFS_DIR)/bin
	chmod +x $(GGNFS_DIR)/bin/gnfs-lasieve4I1*
	touch $(GGNFS_DIR)/bin/gnfs-lasieve4I1*

ggnfs: populate_prefix lasieve_binaries
	cd $(GGNFS_DIR) && \
	$(MAKE) $(MAKE_OPTS) -C $(GGNFS_DIR) x86_64 $(GGNFS_OPTS)

yafu: ecm msieve populate_prefix
	cd $(YAFU_DIR) && \
	$(MAKE) $(MAKE_OPTS) -C $(YAFU_DIR) x86_64 $(YAFU_OPTS)
	cp $(YAFU_DIR)/yafu $(PREFIX)/bin/ 
	cp $(YAFU_DIR)/yafu.ini $(PREFIX)/share/

# All

all: clean populate_prefix ggnfs yafu
	