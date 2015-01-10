YAFU Build Automation Package
=============================

Automate the (grueling) build process for YAFU with NFS support!
This is meant to be used on 64bit systems only. It is likely the case that either the build
or the YAFU binary produced by this will fail on a 32-bit system.

Read more at http://bowser.pw/misc/yafu_setup_package


Building
========

I have written a make script to do all the work for you. Be sure to read caveats at the bottom.

You may specify compile options for the various programs and libraries.

Of interest are 'CUDA=1' and 'OPENMP=1', which may be specified in 'MSIEVE_OPTS',
and 'USE_SSE41=1' which may be specified in 'YAFU_OPTS'.

Once you are ready to build, run `make all`. I would highly suggest that you specify
the maximum threadcount supported by your system, (e.g. `-j4`).

Once compilation is finished, YAFU will be usable inside of the `yafu-XXXXXX` folder, or
can be had in `/prefix/bin`. If you opt to use the YAFU in `/prefix` you will need to 
fill out the `yafu.ini` file properly. You may refer to the `yafu.ini` inside `/prefix/share`
for the correct paths and what have you.


Caveats
=======

I have done my best to put together the makefile and modify the other makefiles as needed.
This compiles without issue on my system, and should take little to no effort on your part.

That being said, please do not attempt to use this from a path that contains spaces.
I looked in to this edge case, but found that mitigating it simply made the other 
Makefiles ***very*** upset.


Happy Factoring, 
	King Bowser.