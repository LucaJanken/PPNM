Out.txt : hello		 	# Out.txt depends on hello
	./hello > Out.txt	# run hello, send output to Out.txt

hello : hello.cc		# hello depends on hello.cc
	g++ -o hello hello.cc	# compile hello.cc into hello

clean: 				# a phoney target, no dependencies
	rm -f Out.txt hello 	# remove secondary files