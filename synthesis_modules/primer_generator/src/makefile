all: primers_gen primers_chk primers_viz

primers_gen: primers_gen.c delta.h display.h filters.h primers.h parameter.h tm.h
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $< -lm
	mv primers_gen ../bin

primers_chk: primers_chk.c primers.h parameter.h filters.h delta.h tm.h index.h hybrid.h display.h
	$(CC) $(CFLAGS) -o $@  $(LDFLAGS) $< -lm
	mv primers_chk ../bin

primers_viz: primers_viz.c primers.h parameter.h filters.h delta.h tm.h index.h hybrid.h display.h
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $< -lm
	mv primers_viz ../bin

