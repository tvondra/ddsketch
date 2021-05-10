MODULE_big = ddsketch
OBJS = ddsketch.o

EXTENSION = ddsketch
DATA = ddsketch--1.0.0.sql
MODULES = ddsketch

CFLAGS=`pg_config --includedir-server`

TESTS        = $(wildcard test/sql/*.sql)
REGRESS      = $(patsubst test/sql/%.sql,%,$(TESTS))
REGRESS_OPTS = --inputdir=test

PG_CONFIG = pg_config
PGXS := $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)
