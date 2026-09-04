# AGENTS.md
Kmers.jl provides the `Kmer` and `Oligo` types for efficient k-mer representation in bioinformatics. This is a BioJulia package that is tightly coupled to BioSequences.jl and relies heavily on its internals.

## Types
* `Kmer{A, K, N}` is an immutable bitstype with `Alphabet` `A`, length K backed by an `NTuple{N, UInt}`. `N` is derived from A and K and is not a free parameter.
  Individual biosymbols (e.g. DNA monomers) are encoded to `UInt`, and then packed from first to last integer of the tuple. Within an integer, they use the lower bits, and are packed from upper to higher bits.
  E.g. a sequence of 16-bit elements "ABCDEFG" would pack as `( ABC, DEFG)`, with the top 16 bits of the first integer being zerod out.

* `Oligo{A, U}` is an immutable bitstype with `Alphabet` `A`, stored in a single unsigned integer of type `U`. The integer encodes its sequence and its length, which is dynamic at runtime. Used to avoid type-instability compared to using `Kmer` when operating on small sequences of varying length.
  Lower bits of the `U` stores the length, the upper bits from top to bottom stores the symbols. Example: `TG` in a `UInt8` would be stored as `0b11_10_00_10`, i.e. T, then G, then unused bits, then the length.

## Coding guidelines
* Performance is paramount in this package, because people use kmers precisely when performance prohibits them from using the more versatile heap-allocated `BioSequences.LongSequence`
* Some functions are intentionally type unstable, such as slicing a `Kmer`, since the length is part of the kmer type. However, this may be type stable if the slice is known at compile time (such as a literal `kmer[1:3]`). For performance reasons, be clear about what code is intentionally type unstable and keep all the rest type stable

* This package relies heavily on zero-cost or cheap abstractions, which rely on inlining for efficiency. The sheer number of helper functions sometimes cause the compiler to not inline, hence the many `@inline` annotations. When creating new functionality that might be performance critical, inspect its generated assembly code and determine if `@inline` is needed.

## Development guidelines
* Format with `runic -i .`. If Runic cannot be run, alert the user instead of finding workarounds for how to run Runic.

* This package depends heavily on BioSequences.jl internals. When in doubt about encoding, check BioSequences code.

* Kmers are optimized for small k. For k > ~100, BioSequences.LongSequence may be more appropriate.

* This package may assume a little-endian 64-bit CPU, and refuse to compile on other platforms.

### RecodingScheme Dispatch
Construction of both Kmer and Oligo uses a `RecodingScheme` abstraction to handle different input types efficiently:

- **`Copyable`**: Direct copy of encoding from compatible BioSequences (same or compatible alphabet)
- **`AsciiEncode`**: Efficient ASCII string parsing for DNA/RNA/AA sequences
- **`TwoToFour`**: Convert 2-bit to 4-bit nucleic acid alphabets
- **`FourToTwo`**: Convert 4-bit to 2-bit nucleic acid alphabets (errors on ambiguous bases)- **`GenericRecoding`**: Fallback for arbitrary iterables

### Alphabet Compatibility
DNA and RNA kmers with the same bit-width (both 2-bit or both 4-bit) are **compatible**:
- They hash to the same value
- They can be compared with `==` and `isequal`
- Construction can convert between them efficiently
