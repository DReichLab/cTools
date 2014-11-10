HOMEL=$(PWD)
TOP=./bin
## binary for install  
# DEBUG_OPTIONS= -g -pg
DEBUG_OPTIONS= -g -pg
BIN=$(HOMEL)/bin
PERLDIR=$(HOMEL)/perlsrc

NLIB=$(HOMEL)/lib/nicklib.a
BLIB=$(HOMEL)/lib/bamlib.a
IDIR=$(HOMEL)/include
##VPATH=.:nicksrc:bamsrc
DEBUG_OPTIONS= -g -pg

CFLAGS= -g -pg -DHAVE_CONFIG_H -I. -I$(IDIR)  -Wimplicit-int \
  -L. -L -lm  -g -pg -D_IOLIB=2  


OBJ=strsubs.o sortit.o vsubs.o statsubs.o linsubs.o getpars.o xsearch.o gauss.o	gds.o
TWTAB=\"$(HOMEL)/tables/twtable\"

##statsubs.o:     tables	statsubs.c
##	$(CC)  $(CFLAGS) -DTWTAB=$(TWTAB) -o statsubs.o statsubs.c

sam_OBJ = bam.o bam_aux.o bam_import.o bam_index.o \
				   bam_lpileup.o bam_maqcns.o bam_mate.o bam_md.o bam_pileup.o \
				   bam_plcmd.o bam_rmdup.o bam_sort.o bam_stat.o bam_rmdupse.o \
				   bam_tview.o bgzf.o faidx.o glf.o kstring.o razf.o sam.o sam_view.o 


IDIR=$(HOMEL)/include  

PROGS=cpulldown cpoly cascertain ccompress cuncompress

Q1I=cpulldown 
Q1IO=cpulldown.o mcio.o admutils.o  $(sam_OBJ) $(OBJ)

allall:all

$(Q1I): nicklib $(Q1IO) 
	gcc -o $(Q1I) $(Q1IO)  -lm -lz -lncurses 

Q2I=cpoly 
Q2IO=cpoly.o $(sam_OBJ) $(OBJ)

$(Q2I): nicklib $(Q2IO) 
	gcc -o $(Q2I) $(Q2IO)  -lm -lz -lncurses 

Q3I=cascertain 
Q3IO=cascertain.o $(sam_OBJ) $(OBJ)

$(Q3I): nicklib $(Q3IO) 
	gcc -o $(Q3I) $(Q3IO)  -lm -lz -lncurses 

Q4I=ccompress 
Q4IO=ccompress.o $(sam_OBJ) $(OBJ)

$(Q4I): nicklib $(Q4IO) 
	gcc -o $(Q4I) $(Q4IO)  -lm -lz -lncurses 

Q5I=cuncompress 
Q5IO=cuncompress.o $(sam_OBJ) $(OBJ)

$(Q5I): nicklib $(Q5IO) 
	gcc -o $(Q5I) $(Q5IO)  -lm -lz -lncurses 

all: $(PROGS) 

install:	all	
	mkdir -p $(TOP)
	cp $(PROGS) $(TOP)

uninstall:
	rm -rf $(NLIB) 
	rm -rf $(TOP)

libnick.a:	dirs $(OBJ)
	ar -r libnick.a $(OBJ)
	ranlib libnick.a

nicklib:	dirs  tables  libnick.a 
	cp libnick.a  $(NLIB)

tables:    
	echo "tables made"  > tables
	cp twtable  $(HOMEL)/tables
	
dirs:	
	mkdir -p  $(HOMEL)/lib
	mkdir -p  $(HOMEL)/include
	mkdir -p  $(BIN)
	cp  *.h  $(IDIR)

clean: 
	rm -f *.o 
	rm -f core
	rm -f libnick.a
	rm -f $(PROGS)

clobber: clean rmdirs uninstall


rmdirs: 
	rm -rf $(HOMEL)/lib 
	rm -rf $(HOMEL)/include 

