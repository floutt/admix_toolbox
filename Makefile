CC       ?= cc

DOTGENO_PREFIX ?= /usr/local

PREFIX ?= /usr/local
BINDIR  = $(PREFIX)/bin

SRCDIR = src
OUTDIR = bin

CPPFLAGS += -I$(DOTGENO_PREFIX)/include -D_POSIX_C_SOURCE=200809L
CFLAGS   ?= -std=c17 -Wall -Wextra -Wpedantic -O2
LDFLAGS  += -L$(DOTGENO_PREFIX)/lib -Wl,-rpath,$(DOTGENO_PREFIX)/lib
LDLIBS   += -ldotgeno -lm

SRCS  = $(wildcard $(SRCDIR)/*.c)
BINS  = $(patsubst $(SRCDIR)/%.c,$(OUTDIR)/%,$(SRCS))

.PHONY: all clean install uninstall

all: $(BINS)

$(OUTDIR):
	mkdir -p $(OUTDIR)

$(OUTDIR)/%: $(SRCDIR)/%.c | $(OUTDIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -o $@ $< $(LDLIBS)

install: all
	mkdir -p $(DESTDIR)$(BINDIR)
	install -m 755 $(BINS) $(DESTDIR)$(BINDIR)/

uninstall:
	rm -f $(addprefix $(DESTDIR)$(BINDIR)/,$(notdir $(BINS)))

clean:
	rm -rf $(OUTDIR)
