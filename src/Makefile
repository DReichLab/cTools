HOMEL=$(PWD)
TOP=./bin
## binary for install  
DEBUG_OPTIONS= -g #-pg
BIN=$(HOMEL)/bin
PERLDIR=$(HOMEL)/perlsrc

NLIB=$(HOMEL)/lib/nicklib.a
BLIB=$(HOMEL)/lib/bamlib.a
IDIR=$(HOMEL)/include
##VPATH=.:nicksrc:bamsrc

CFLAGS= -g  -DHAVE_CONFIG_H -I. -I$(IDIR)  -Wimplicit-int \
   -L. -lm  -g -D_IOLIB=2 #-O3 


OBJ= nicklib/strsubs.o nicklib/sortit.o nicklib/vsubs.o nicklib/statsubs.o nicklib/linsubs.o nicklib/getpars.o nicklib/xsearch.o nicklib/gauss.o	nicklib/gds.o
TWTAB=\"$(HOMEL)/tables/twtable\"

##statsubs.o:     tables	statsubs.c
##	$(CC)  $(CFLAGS) -DTWTAB=$(TWTAB) -o statsubs.o statsubs.c

sam_OBJ = samtools_api/bam.o samtools_api/bam_aux.o samtools_api/bam_import.o samtools_api/bam_index.o \
				   samtools_api/bam_lpileup.o samtools_api/bam_maqcns.o samtools_api/bam_mate.o samtools_api/bam_md.o samtools_api/bam_pileup.o \
				   samtools_api/bam_plcmd.o samtools_api/bam_rmdup.o samtools_api/bam_sort.o samtools_api/bam_stat.o samtools_api/bam_rmdupse.o \
				   samtools_api/bam_tview.o samtools_api/bgzf.o samtools_api/faidx.o samtools_api/glf.o samtools_api/kstring.o samtools_api/razf.o samtools_api/sam.o samtools_api/sam_view.o 


IDIR=$(HOMEL)/include  

PROGS=cpulldown cpoly cascertain ccompress puncompress cmakefilter filtstats

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

Q5I=puncompress 
Q5IO=uncomp.o $(sam_OBJ) $(OBJ)

$(Q5I): nicklib $(Q5IO) 
	gcc -o $(Q5I) $(Q5IO)  -lm -lz -lncurses 

Q6I=cmakefilter
Q6IO=cmakefilter.o opth.o $(sam_OBJ) $(OBJ)

$(Q6I): nicklib $(Q6IO)
	gcc -o $(Q6I) $(Q6IO) -lm -lz -lncurses
#    gcc -o $(Q6I) $(Q6IO)  -lm -lz -lncurses

Q7I=filtstats
Q7IO=filtstats.o $(sam_OBJ) $(OBJ)

$(Q7I): nicklib $(Q7IO)
	gcc -o $(Q7I) $(Q7IO) -lm -lz -lncurses
#    gcc -o $(Q7I) $(Q7IO)  -lm -lz -lncurses

all: $(PROGS) 

install:	all	
	mkdir -p $(TOP)
	cp $(PROGS) uncompress.pl $(TOP)

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
	rm -f samtools_api/*.o 
	rm -f nicklib/*.o 
	rm -f core
	rm -f libnick.a
	rm -f $(PROGS)

clobber: clean rmdirs uninstall


rmdirs: 
	rm -rf $(HOMEL)/lib 
	rm -rf $(HOMEL)/include 

