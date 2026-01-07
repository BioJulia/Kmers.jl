# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

Kmers.jl provides the `Kmer` and `DynamicKmer` types for efficient k-mer representation in bioinformatics. This is a BioJulia package that is tightly coupled to BioSequences.jl and relies heavily on its internals.

**Key Concepts:**
- **Kmer**: Immutable bitstype sequences of fixed length `K` (compile-time). Stored directly in registers for maximum efficiency. Parameterized as `Kmer{A,K,N}` where `A` is the Alphabet, `K` is the length, and `N` is the number of UInt tuples (derived, not a free parameter).
- **DynamicKmer**: Similar to Kmer but with runtime length. Slightly less efficient but useful when working with varying k-mer sizes to avoid excessive compilation and type instability.

**Performance Warning:** Kmers are highly optimized for small, fixed-length sequences. Operations that change length (push, pop, slicing) are type unstable unless the compiler can use constant folding. Kmers become inefficient for longer sequences (e.g., reverse-complement of 512-mer takes 16 μs vs 126 ns for LongSequence).

## Development Commands

### Running Tests

```bash
# Run all tests
julia --project -e 'using Pkg; Pkg.test()'

# Run tests from REPL
julia --project
using Pkg
Pkg.test()

# Run specific test file
julia --project test/dynamic.jl

# Run tests with coverage (if needed)
julia --project --code-coverage test/runtests.jl
```

### Building Documentation

```bash
cd docs/
julia --project -e 'using Pkg; Pkg.instantiate(); include("make.jl")'
```

### Julia Package Development

```bash
# Activate the project
julia --project

# Install dependencies
using Pkg
Pkg.instantiate()

# Run the package
using Kmers
```

## Architecture

### Core Type System

**Kmer Layout (`struct Kmer{A,K,N}`):**
- Data stored as `NTuple{N, UInt}`
- Symbols pack into the **lowest bits** of each UInt, with **first symbols in higher parts**
- Example: A 16-bit element sequence "A-G" would pack as `(ABC, DEFG)` with 16 unused bits at top
- Unused bits are always **zero** and always in the **top bits of first UInt**
- This layout simplifies comparison operators at the cost of more complex construction

**DynamicKmer Layout (`struct DynamicKmer{A,U}`):**
- Single unsigned integer `x::U` containing both data and length
- Lower bits store the **length**
- Upper bits store the **sequence data** (from top to bottom)
- Example with 2-bit alphabet and UInt8: `TG` is stored as `11 10 00 10` (T G unused length=2)

### RecodingScheme Dispatch

Construction of both Kmer and DynamicKmer uses a `RecodingScheme` pattern to handle different input types efficiently:

- **`Copyable`**: Direct copy from compatible BioSequences (same or compatible alphabet)
- **`AsciiEncode`**: Efficient ASCII string parsing for DNA/RNA/AA sequences
- **`TwoToFour`**: Convert 2-bit to 4-bit nucleic acid alphabets
- **`FourToTwo`**: Convert 4-bit to 2-bit nucleic acid alphabets (errors on ambiguous bases)
- **`GenericRecoding`**: Fallback for arbitrary iterables

### Source Structure

- **`src/kmer.jl`**: Core Kmer type definition, type checking, basic operations
- **`src/dynamic.jl`**: DynamicKmer type with runtime length
- **`src/construction.jl`**: Kmer construction logic and RecodingScheme dispatch
- **`src/construction_utils.jl`**: Unsafe extraction and shifting utilities for iterators
- **`src/indexing.jl`**: Indexing operations (scalar, range, logical)
- **`src/transformations.jl`**: Biological operations (reverse, complement, canonical)
- **`src/revtrans.jl`**: Reverse translation (amino acid → codons)
- **`src/counting.jl`**: GC content, symbol counting
- **`src/tuple_bitflipping.jl`**: Low-level bit manipulation for NTuple operations
- **`src/iterators/`**: Various kmer iterators (FwKmers, CanonicalKmers, UnambiguousKmers, SpacedKmers)

### Test Structure

Tests are organized by feature:
- **`test/runtests.jl`**: Main test runner with comprehensive Kmer tests
- **`test/dynamic.jl`**: DynamicKmer-specific tests (included from runtests.jl)
- **`test/translation.jl`**: Translation and genetic code tests
- **`test/benchmark.jl`**: Performance benchmarks
- **`test/utils.jl`**: Test utilities

## Important Implementation Details

### Type Parameters and Derivation

The `N` parameter in `Kmer{A,K,N}` is **not a free parameter** - it's derived from `A` and `K` using `derive_type`:
```julia
derive_type(Kmer{DNAAlphabet{2}, 5})  # Returns Kmer{DNAAlphabet{2}, 5, 1}
```

When constructing kmers or writing tests, you can omit `N` and it will be derived automatically.

### Alphabet Compatibility

DNA and RNA kmers with the same bit-width (both 2-bit or both 4-bit) are **compatible**:
- They hash to the same value
- They can be compared with `==` and `isequal`
- Construction can convert between them efficiently

However, they are **distinct types** and cannot be directly compared with other BioSequences (this throws an exception).

### Integer Conversion

- `as_integer(kmer)`: Extract coding bits only (no length information)
- `from_integer(Type, u)`: Reconstruct from integer (for Kmer)
- `from_integer(Type, u, len)`: Reconstruct from integer with explicit length (for DynamicKmer)

For DynamicKmer, `as_integer` and `from_integer` with different lengths may produce reproducible but incorrect results.

### Unsafe Methods

Several internal functions use the `Unsafe` trait object (`unsafe` singleton) to bypass bounds checking:
- `unsafe_extract`: Extract kmer from sequence without bounds checking
- `unsafe_shift_from`: Shift symbols into kmer from sequence
- These are used internally by iterators for performance

### Writing Tests for New Features

When adding tests, follow the existing patterns:
1. Group related tests in `@testset` blocks with descriptive names
2. Test edge cases: empty sequences, length=1, maximum length
3. Test multiple alphabets where applicable (DNA 2-bit, DNA 4-bit, RNA, AA)
4. Test different integer widths for DynamicKmer (UInt32, UInt64, UInt128)
5. Test round-trip conversions (construct → as_integer → from_integer)
6. For DynamicKmer: test widening and narrowing of backing integer types

## Common Patterns

### Testing Kmer Construction from Various Sources

```julia
for s in [dna"TAGCTA", rna"UGCUGA", aa"PLKWM"]
    kmer = Kmer{typeof(Alphabet(s)), length(s)}(s)
    @test string(kmer) == string(s)
end
```

### Testing DynamicKmer Type Conversions

```julia
# Same alphabet, different backing type
d32 = DynamicDNAKmer{UInt32}(dna"TAGC")
d64 = DynamicDNAKmer{UInt64}(d32)
@test d64 == d32

# Different alphabets (DNA/RNA)
d_dna = DynamicDNAKmer{UInt64}(dna"ATGT")
d_rna = DynamicRNAKmer{UInt64}(d_dna)
@test d_rna == d_dna
```

## Notes for AI Assistants

1. **BioSequences Integration**: This package depends heavily on BioSequences.jl internals. When in doubt about encoding, check BioSequences documentation.

2. **Type Stability**: Length-changing operations on Kmer are inherently type-unstable. Consider recommending DynamicKmer for variable-length use cases.

3. **Testing Philosophy**: The test suite is comprehensive. New features should have similar test coverage including edge cases, multiple alphabets, and round-trip conversions.

4. **Performance**: Kmers are optimized for small k. For k > ~100, LongSequence may be more appropriate.

5. **Little-Endian Assumption**: Some code (particularly Kmer to DynamicKmer conversion) assumes little-endian architecture. The package won't compile on big-endian systems.
