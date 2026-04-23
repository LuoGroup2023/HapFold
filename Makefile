SRC_DIR = src
LIB_DIR = lib

CFLAGS = -g -gdwarf-3 -fpermissive -Wall -O0 #-O2
CC = g++


INCLUDES = -I. -I$(SRC_DIR)


UTILS_OBJS = hash.o dict.o array.o utils.o
CORE_OBJS = kthread.o bbf.o htab.o bseq.o misc.o kalloc.o paf.o \
            hic_mapping.o hic_completeness.o hic_qv.o count.o \
            hic_switch_error.o seqio.o seqhash.o $(UTILS_OBJS)


OBJS = $(addprefix $(SRC_DIR)/, $(CORE_OBJS))


PROG = HapFold


LIBS = -lm -lz -lpthread $(LIB_DIR)/libminimap2.a $(LIB_DIR)/libz.a

ifneq ($(asan),)
    CFLAGS += -fsanitize=address
    LIBS += -fsanitize=address
endif

.PHONY: all clean

all: $(PROG)


$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@


$(SRC_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@


main.o: main.cpp
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@


$(PROG): $(OBJS) main.o
	$(CC) $(CFLAGS) $^ $(LIBS) -o $@

clean:
	rm -fr gmon.out $(SRC_DIR)/*.o *.o ext/*.o a.out $(PROG) *~ *.dSYM session*