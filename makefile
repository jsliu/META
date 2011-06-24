CFLAGS = -O3 -static  
CXXFLAGS = -I. -O3 -Wno-deprecated -static

CXX = /usr/bin/g++
AR = /usr/bin/ar

PROB_LIBS = -L./tools/cprob/ -I./tools/cprob/ 
CPROB := $(patsubst %.c,%.o,$(wildcard tools/cprob/*.c))

#source files
SOURCE = main.o meta_analysis.o

all: meta
meta: $(SOURCE) libcprob
	$(CXX) $(CXXFLAGS) $(SOURCE) -o meta-stat-1.3 -lz -lm -lcprob $(PROB_LIBS) -static

$(SOURCE): meta_analysis.h libcprob
$(BOOST_IO): %.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< 
libcprob: $(CPROB)
	$(AR) rc tools/cprob/libcprob.a $(CPROB)

clean:
	rm -f *.o ./tools/cprob/*.o ./tools/cprob/*.a; rm meta-*
