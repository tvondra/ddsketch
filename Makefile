# The version the distribution archive and the release notes are built for.
EXTENSION    = $(shell grep -m 1 '"name":' META.json | \
               sed -e 's/[[:space:]]*"name":[[:space:]]*"\([^"]*\)",/\1/')
EXTVERSION   = $(shell grep -m 1 '[[:space:]]\{8\}"version":' META.json | \
               sed -e 's/[[:space:]]*"version":[[:space:]]*"\([^"]*\)",\{0,1\}/\1/')
DISTVERSION  = $(shell grep -m 1 '[[:space:]]\{3\}"version":' META.json | \
               sed -e 's/[[:space:]]*"version":[[:space:]]*"\([^"]*\)",\{0,1\}/\1/')

MODULE_big = ddsketch
OBJS = ddsketch.o

DATA = sql/ddsketch--1.0.0.sql sql/ddsketch--1.0.0--1.0.1.sql sql/ddsketch--1.0.1--1.0.2-dev.sql sql/ddsketch--1.1.0-dev.sql

TESTS        = $(wildcard test/sql/*.sql)
REGRESS      = $(patsubst test/sql/%.sql,%,$(TESTS))
REGRESS_OPTS = --inputdir=test

PG_CONFIG = pg_config
PGXS := $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)

# The PGXS "uninstall" target only removes the files currently listed in DATA,
# so scripts installed by an older build (e.g. upgrade scripts that were since
# renamed or removed) would be left behind. Match the control file and versioned
# SQL scripts, not just the name prefix, which other extensions may share.
#
# A recipe can't be appended to the PGXS "uninstall" rule, so this is hooked in
# as a prerequisite (which means it runs before the PGXS part - that's fine,
# the two steps are independent).
uninstall: uninstall-stale

.PHONY: uninstall-stale
uninstall-stale:
	rm -f '$(DESTDIR)$(datadir)/extension/$(EXTENSION).control' \
		'$(DESTDIR)$(datadir)/extension/$(EXTENSION)--'*.sql

dist:
	git archive --format zip --prefix=$(EXTENSION)-$(DISTVERSION)/ -o $(EXTENSION)-$(DISTVERSION).zip HEAD

# The output is tracked, so checkout timestamps do not establish whether its
# contents match the release selected by META.json.
.PHONY: latest-changes.md
latest-changes.md: Changes META.json
	perl -e 'while (<>) {last if /^(v?\Q${DISTVERSION}\E)/; } print "Changes for v${DISTVERSION}:\n"; while (<>) { last if /^\s*$$/; s/^\s+//; print }' Changes > $@
