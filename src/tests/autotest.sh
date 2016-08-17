# /opt/local/include points to boost
./configure CPPFLAGS=-I/opt/local/include && make check
./configure CXX=clang++-mp-3.7  CPPFLAGS=-I/opt/local/include  && make check
