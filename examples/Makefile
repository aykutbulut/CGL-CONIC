CC=g++
PROGNAME=conic_cgl
INCLUDES=`pkg-config --cflags cola cglconic`
CFLAGS=-g
#LDFLAGS=-L/scratch/aykut/alps/build/lib -L/usr/lib/ -L/usr/local/coinall/1.7.4/lib
LIBS=`pkg-config --libs cola cglconic`

OBJECTS=conic_cgl.o

all: $(PROGNAME)

$(PROGNAME): $(OBJECTS)
	$(CC) -o $(PROGNAME) $(OBJECTS)  $(LDFLAGS) $(LIBS)

$(OBJECTS): Makefile

.cpp.o:
	$(CC) -c $(CFLAGS) $(INCLUDES) -o $@ $<

clean:
	rm *.o $(PROGNAME)

test: $(PROGNAME)
	./$(PROGNAME) ~/cblib/estein4_A.mps min
