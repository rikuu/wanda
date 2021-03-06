# Wanda

Wanda is a fast, low-memory genome assembler toy.

## Usage

```sh
$ concatenate <output> <file> # concatenates sequences into a stream-like format
$ wanda-build <stream> <k> <graph prefix> # builds indices
$ wanda-assemble <graph prefix> <s> <min length> # assembles unitigs
```

## Dependencies
- A compiler that supports C++11,
- [SDSL-lite][sdsl-lite] - low level succinct data structures,
- [pSAscan][psascan] - parallel external memory suffix array construction.

[sdsl-lite]: https://github.com/simongog/sdsl-lite
[psascan]: https://www.cs.helsinki.fi/group/pads/pSAscan.html
