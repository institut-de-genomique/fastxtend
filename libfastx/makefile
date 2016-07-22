all : libfastx.a

chomp.o : .deps/chomp.Po chomp.c chomp.h
	if gcc -DHAVE_CONFIG_H -I. -I. -I../..     -g -O2 -Wall -Wextra -Wformat-nonliteral -Wformat-security -Wswitch-default -Wswitch-enum -Wunused-parameter -Wfloat-equal -DDEBUG -g -O1 -MT chomp.o -MD -MP -MF ".deps/chomp.Tpo" -c -o chomp.o chomp.c; \
	then mv -f ".deps/chomp.Tpo" ".deps/chomp.Po"; else rm -f ".deps/chomp.Tpo"; exit 1; fi

fastx.o : .deps/fastx.Po fastx.c fastx.h
	if gcc -DHAVE_CONFIG_H -I. -I. -I../..     -g -O2 -Wall -Wextra -Wformat-nonliteral -Wformat-security -Wswitch-default -Wswitch-enum -Wunused-parameter -Wfloat-equal -DDEBUG -g -O1 -MT fastx.o -MD -MP -MF ".deps/fastx.Tpo" -c -o fastx.o fastx.c; \
	then mv -f ".deps/fastx.Tpo" ".deps/fastx.Po"; else rm -f ".deps/fastx.Tpo"; exit 1; fi

fastx_args.o : .deps/fastx_args.Po fastx_args.c fastx_args.h
	if gcc -DHAVE_CONFIG_H -I. -I. -I../..     -g -O2 -Wall -Wextra -Wformat-nonliteral -Wformat-security -Wswitch-default -Wswitch-enum -Wunused-parameter -Wfloat-equal -DDEBUG -g -O1 -MT fastx_args.o -MD -MP -MF ".deps/fastx_args.Tpo" -c -o fastx_args.o fastx_args.c; \
	then mv -f ".deps/fastx_args.Tpo" ".deps/fastx_args.Po"; else rm -f ".deps/fastx_args.Tpo"; exit 1; fi

sequence_alignment.o : .deps/sequence_alignment.Po sequence_alignment.cpp sequence_alignment.h
	if g++ -DHAVE_CONFIG_H -I. -I. -I../..     -g -O2 -Wall -Wextra -Wformat-nonliteral -Wformat-security -Wswitch-default -Wswitch-enum -Wunused-parameter -Wfloat-equal -DDEBUG -g -O1 -DDEBUG -g -O1 -MT sequence_alignment.o -MD -MP -MF ".deps/sequence_alignment.Tpo" -c -o sequence_alignment.o sequence_alignment.cpp; \
	then mv -f ".deps/sequence_alignment.Tpo" ".deps/sequence_alignment.Po"; else rm -f ".deps/sequence_alignment.Tpo"; exit 1; fi

libfastx.a : chomp.o fastx.o fastx_args.o sequence_alignment.o
	rm -f libfastx.a; ar cru libfastx.a chomp.o fastx.o fastx_args.o sequence_alignment.o; ranlib libfastx.a;

clean :
	rm -f libfastx.a chomp.o fastx.o fastx_args.o sequence_alignment.o
