# SRC_DIR = src
# LIB_DIR = lib

# CFLAGS = -g -gdwarf-3 -fpermissive -Wall -O0 #-O2
# CC = g++


# INCLUDES = -I. -I$(SRC_DIR)


# UTILS_OBJS = hash.o dict.o array.o utils.o
# CORE_OBJS = kthread.o bbf.o htab.o bseq.o misc.o kalloc.o paf.o \
#             mapping.o seqio.o seqhash.o $(UTILS_OBJS)


# OBJS = $(addprefix $(SRC_DIR)/, $(CORE_OBJS))


# PROG = HapFold


# LIBS = -lm -lz -lpthread $(LIB_DIR)/libminimap2.a $(LIB_DIR)/libz.a

# ifneq ($(asan),)
#     CFLAGS += -fsanitize=address
#     LIBS += -fsanitize=address
# endif

# .PHONY: all clean

# all: $(PROG)


# $(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
# 	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@


# $(SRC_DIR)/%.o: $(SRC_DIR)/%.c
# 	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@


# main.o: main.cpp
# 	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@


# $(PROG): $(OBJS) main.o
# 	$(CC) $(CFLAGS) $^ $(LIBS) -o $@

# clean:
# 	rm -fr gmon.out $(SRC_DIR)/*.o *.o ext/*.o a.out $(PROG) *~ *.dSYM session*



# ============================================================
# Project directories
# ============================================================

SRC_DIR       := src
HIFIASM_DIR   := hifiasm
LIB_DIR       := lib
BUILD_DIR     := build

HF_OBJ_DIR    := $(BUILD_DIR)/hapfold
HA_OBJ_DIR    := $(BUILD_DIR)/hifiasm

PROG          := HapFold


# ============================================================
# Compilers and tools
# ============================================================

CXX           := g++
CC            := gcc
AR            := ar

#
# Preserve HapFold's existing use of g++ for .c files to avoid changing
# established C/C++ linkage behavior.
#
HF_CC         := $(CXX)


# ============================================================
# Compile flags
# ============================================================

HF_CXXFLAGS   := -O3 -std=c++17 -fpermissive -Wall
HA_CXXFLAGS   := -g -O3 -std=c++17 -msse4.2 -mpopcnt \
                 -fomit-frame-pointer -Wall -fPIC -fvisibility=hidden

DEPFLAGS      := -MMD -MP

HF_INCLUDES   := -I. -I$(SRC_DIR)
HA_INCLUDES   := -I$(HIFIASM_DIR)


# ============================================================
# HapFold objects
# ============================================================

HF_UTIL_NAMES := hash.o dict.o array.o utils.o

HF_CORE_NAMES := \
	kthread.o \
	bbf.o \
	htab.o \
	bseq.o \
	misc.o \
	kalloc.o \
	paf.o \
	mapping.o \
	seqio.o \
	seqhash.o \
	sys.o \
	$(HF_UTIL_NAMES)

HF_OBJS       := $(addprefix $(HF_OBJ_DIR)/,$(HF_CORE_NAMES))
HF_MAIN_OBJ   := $(HF_OBJ_DIR)/main.o


# ============================================================
# hifiasm objects
# ============================================================

HA_ALL_NAMES := \
	CommandLines.o \
	Process_Read.o \
	Assembly.o \
	Hash_Table.o \
	POA.o \
	Correct.o \
	Levenshtein_distance.o \
	Overlaps.o \
	Trio.o \
	kthread.o \
	Purge_Dups.o \
	htab.o \
	hist.o \
	sketch.o \
	anchor.o \
	extract.o \
	sys.o \
	hic.o \
	rcut.o \
	horder.o \
	ecovlp.o \
	tovlp.o \
	inter.o \
	kalloc.o \
	gfa_ut.o \
	gchain_map.o

#
# Hifiasm uses its complete native implementations. Same-named source files
# are not assumed to be ABI-compatible, so a hidden-symbol shared library
# isolates hifiasm globals from HapFold globals.
#
HA_OBJS       := $(addprefix $(HA_OBJ_DIR)/,$(HA_ALL_NAMES))
HA_MAIN_OBJ   := $(HA_OBJ_DIR)/main.o

HA_SHARED_LIB := $(BUILD_DIR)/libhifiasm_embedded.so


# ============================================================
# Libraries
# ============================================================

#
# Avoid linking both -lz and lib/libz.a, which can mix zlib versions.
# To use the system zlib instead, replace the definition below with:
#
# LIBS := $(LIB_DIR)/libminimap2.a -lz -lpthread -lm
#
LIBS := \
	$(LIB_DIR)/libminimap2.a \
	$(LIB_DIR)/libz.a \
	-lpthread \
	-lm


# ============================================================
# AddressSanitizer
# ============================================================

ifneq ($(asan),)
	HF_CXXFLAGS += -fsanitize=address \
	               -fno-omit-frame-pointer
	HA_CXXFLAGS += -fsanitize=address \
	               -fno-omit-frame-pointer
	LIBS        += -fsanitize=address
endif


# ============================================================
# Targets
# ============================================================

.PHONY: all clean check-shared standalone-hifiasm

all: $(PROG)


# ============================================================
# Build directories
# ============================================================

$(BUILD_DIR):
	mkdir -p $@

$(HF_OBJ_DIR):
	mkdir -p $@

$(HA_OBJ_DIR):
	mkdir -p $@


# ============================================================
# HapFold compilation rules
# ============================================================

$(HF_OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(HF_OBJ_DIR)
	$(HF_CC) \
		$(HF_CXXFLAGS) \
		$(DEPFLAGS) \
		$(HF_INCLUDES) \
		-c $< \
		-o $@

$(HF_OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(HF_OBJ_DIR)
	$(HF_CC) \
		$(HF_CXXFLAGS) \
		$(DEPFLAGS) \
		$(HF_INCLUDES) \
		-c $< \
		-o $@

$(HF_MAIN_OBJ): main.cpp \
                $(HIFIASM_DIR)/hifiasm_entry.h \
                | $(HF_OBJ_DIR)
	$(CXX) \
		$(HF_CXXFLAGS) \
		$(DEPFLAGS) \
		$(HF_INCLUDES) \
		-c $< \
		-o $@


# ============================================================
# hifiasm compilation rules
# ============================================================

$(HA_OBJ_DIR)/%.o: $(HIFIASM_DIR)/%.cpp | $(HA_OBJ_DIR)
	$(CXX) \
		$(HA_CXXFLAGS) \
		$(DEPFLAGS) \
		$(HA_INCLUDES) \
		-c $< \
		-o $@

$(HA_OBJ_DIR)/%.o: $(HIFIASM_DIR)/%.c | $(HA_OBJ_DIR)
	$(CC) \
		$(HA_CXXFLAGS) \
		$(DEPFLAGS) \
		$(HA_INCLUDES) \
		-c $< \
		-o $@

$(HA_MAIN_OBJ): $(HIFIASM_DIR)/hifiasm_main.cpp \
                $(HIFIASM_DIR)/hifiasm_entry.h \
                | $(HA_OBJ_DIR)
	$(CXX) \
		$(HA_CXXFLAGS) \
		$(DEPFLAGS) \
		$(HA_INCLUDES) \
		-c $< \
		-o $@


# ============================================================
# Embedded hifiasm shared library
# ============================================================

$(HA_SHARED_LIB): $(HA_OBJS) $(HA_MAIN_OBJ) | $(BUILD_DIR)
	$(CXX) -shared -Wl,-Bsymbolic \
		$(HA_OBJS) $(HA_MAIN_OBJ) \
		-lz -lpthread -lm \
		-o $@


# ============================================================
# Final HapFold executable
# ============================================================

$(PROG): \
	$(HF_OBJS) \
	$(HF_MAIN_OBJ) \
	$(HA_SHARED_LIB)

	$(CXX) \
		$(HF_CXXFLAGS) \
		$(HF_OBJS) \
		$(HF_MAIN_OBJ) \
		-Wl,--start-group \
		$(LIBS) \
		-Wl,--end-group \
		-L$(BUILD_DIR) -lhifiasm_embedded \
		-Wl,-rpath,'$$ORIGIN/$(BUILD_DIR)' \
		-o $@


# ============================================================
# Compare shared modules before linking
# ============================================================

check-shared:
	@echo "===== kthread headers ====="
	-diff -u $(SRC_DIR)/kthread.h \
	         $(HIFIASM_DIR)/kthread.h

	@echo "===== kalloc headers ====="
	-diff -u $(SRC_DIR)/kalloc.h \
	         $(HIFIASM_DIR)/kalloc.h

	@echo "===== htab headers ====="
	-diff -u $(SRC_DIR)/htab.h \
	         $(HIFIASM_DIR)/htab.h

	@echo "===== source checksums ====="
	-sha256sum \
		$(SRC_DIR)/kthread.c \
		$(HIFIASM_DIR)/kthread.cpp \
		$(SRC_DIR)/kalloc.c \
		$(HIFIASM_DIR)/kalloc.cpp \
		$(SRC_DIR)/htab.cpp \
		$(HIFIASM_DIR)/htab.cpp \
		$(SRC_DIR)/sys.cpp \
		$(HIFIASM_DIR)/sys.cpp


# ============================================================
# Optional standalone hifiasm
# ============================================================

standalone-hifiasm:
	$(MAKE) \
		-C $(HIFIASM_DIR) \
		CXXFLAGS="$(HA_CXXFLAGS) -DHIFIASM_STANDALONE"


# ============================================================
# Dependencies
# ============================================================

DEPS := \
	$(HF_OBJS:.o=.d) \
	$(HF_MAIN_OBJ:.o=.d) \
	$(HA_OBJS:.o=.d) \
	$(HA_MAIN_OBJ:.o=.d)

-include $(DEPS)


# ============================================================
# Clean
# ============================================================

clean:
	rm -rf \
		$(BUILD_DIR) \
		$(PROG) \
		gmon.out \
		a.out \
		*.dSYM \
		session* \
		*~
