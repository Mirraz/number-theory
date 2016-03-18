CC=g++
LD=g++
STRIP=strip -s
WARNINGS=-Wall -Wextra -pedantic
DEBUG=
#DEBUG=-g -ggdb
COPTIM=-march=native -O2
#COPTIM=-O0
DEFINES=
INCLUDES=
CSTD=-std=c++11
CFLAGS=$(WARNINGS) $(DEBUG) $(COPTIM) $(DEFINES) $(INCLUDES) $(CSTD) -pipe
LDOPTIM=-Wl,-O1 -Wl,--as-needed
#LDOPTIM=
LIBFILES=-lm
LDFLAGS=$(WARNINGS) $(DEBUG) $(LDOPTIM) $(LIBFILES)
SRC_DIR=.
BUILD_DIR=build

ALL_TESTS=factorize_tests pow_mod_tests primitive_roots_tests

tests: $(ALL_TESTS)

factorize_tests: $(BUILD_DIR)/factorize_tests.o
	$(LD) -o $@ $^ $(LDFLAGS)
	$(STRIP) $@

$(BUILD_DIR)/factorize_tests.o: $(SRC_DIR)/factorize_tests.cpp $(SRC_DIR)/factorize.h Makefile
	$(CC) -o $@ $< -c $(CFLAGS)

pow_mod_tests: $(BUILD_DIR)/pow_mod_tests.o
	$(LD) -o $@ $^ $(LDFLAGS)
	$(STRIP) $@

$(BUILD_DIR)/pow_mod_tests.o: $(SRC_DIR)/pow_mod_tests.cpp $(SRC_DIR)/pow_mod.h Makefile
	$(CC) -o $@ $< -c $(CFLAGS)

primitive_roots_tests: $(BUILD_DIR)/primitive_roots_tests.o
	$(LD) -o $@ $^ $(LDFLAGS)
	$(STRIP) $@

$(BUILD_DIR)/primitive_roots_tests.o: $(SRC_DIR)/primitive_roots_tests.cpp $(SRC_DIR)/pow_mod.h $(SRC_DIR)/factorize.h Makefile
	$(CC) -o $@ $< -c $(CFLAGS)

clean_tests:
	rm $(ALL_TESTS)

clean:
	rm -r $(BUILD_DIR)
	mkdir $(BUILD_DIR)

