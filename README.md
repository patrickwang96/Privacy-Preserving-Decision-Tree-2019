# Privacy-Preserving-Decision-Tree-2019

## Network Setup

IP address is hard coded into network.h file. Default is localhost. Please change the ip address before compilation.

## Compile

```bash
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../

make server client -j4
```

