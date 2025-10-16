# DNA Extension Makefile

MODULE_big = dna_ext
OBJS = \
	src/module.o \
	src/dna_utils.o \
	src/type_dna.o \
	src/type_kmer.o \
	src/type_qkmer.o \
	src/funcs.o \
	src/ops.o \
	src/hash_ops.o \
	src/btree_ops.o
#	src/spgist_kmer.o

EXTENSION = dna_ext
DATA = sql/dna_ext--1.0.sql
PGFILEDESC = "dna_ext - DNA sequence data types for PostgreSQL"

# Include PostgreSQL extension build system
ifdef USE_PGXS
PG_CONFIG = pg_config
PGXS := $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)
else
subdir = contrib/dna_ext
top_builddir = ../..
include $(top_builddir)/src/Makefile.global
include $(top_srcdir)/contrib/contrib-global.mk
endif

# Additional compilation flags
CFLAGS += -Wall -Wextra -std=c99

# Clean target
clean: 
	rm -f $(OBJS)
	rm -f dna_ext.so

.PHONY: clean install uninstall
