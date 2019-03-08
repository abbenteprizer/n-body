# Usage:
# make        # compile all binary
# make clean  # remove ALL binaries and objects

.PHONY	= all clean

CC = g++ # compiler to use

LINKERFLAG = -lm

SRCS := snbody.cpp
#SRCS := pnbody.cpp
BINS := $(SRCS:%.cpp=%)

all:	${BINS}

%: %.o
	@echo "Checking.."
	${CC}	${LINKERFLAG}	$<-o$@

%.o: %.c
	@echo "Creating object.."
	${CC}	-c	$<

:	./snbody > hello.dat

clean:
	@echo "Cleaning up..."
	rm -rvf *.o	${BINS}

# run and save data

# create image using gnuplot
