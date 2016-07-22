OBJLIBS = ../libfastx/libfastx.a
 
all : fastx_stats

fastx_stats : fastx_stats.cpp .deps/fastx_stats.Po $(OBJLIBS)
	if g++ -DHAVE_CONFIG_H -I. -I. -I../..   -I../libfastx   -g -O2 -Wall -Wextra -Wformat-nonliteral -Wformat-security -Wswitch-default -Wswitch-enum -Wunused-parameter -Wfloat-equal -DDEBUG -g -O1 -DDEBUG -g -O1 -MT fastx_stats.o -MD -MP -MF ".deps/fastx_stats.Tpo" -c -o fastx_stats.o fastx_stats.cpp; \
	then mv -f ".deps/fastx_stats.Tpo" ".deps/fastx_stats.Po"; else rm -f ".deps/fastx_stats.Tpo"; exit 1; fi; \
	g++  -g -O2 -Wall -Wextra -Wformat-nonliteral -Wformat-security -Wswitch-default -Wswitch-enum -Wunused-parameter -Wfloat-equal -DDEBUG -g -O1 -DDEBUG -g -O1   -o fastx_stats  fastx_stats.o ../libfastx/libfastx.a

install :
	cp -f fastx_stats ../bin

clean :
	rm -f fastx_stats fastx_stats.o
