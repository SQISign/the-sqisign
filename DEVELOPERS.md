# Developer guidelines

Please read carefully before contributing to this repo.

## Code structure

The code is split into the modules below:

- `common`: common code for AES, SHAKE, (P)RNG, memory handling. Every
  module that needs a hash function, seed expansion (e.g., KLPT),
  deterministic alea for tests, should call to this module.
- `uintbig`: multi-precision big integers.
- `gf`: GF(p^2) and GF(p) arithmetic.
- `ec`: elliptic curves, isogenies and pairings. Everything that is
  purely finite-fieldy.
- `quaternion`: quaternion orders and ideals. This is, essentially,
  replacing PARI/GP.
- `klpt`: implementation of KLPT.
- `id2iso`: code for Iso <-> Ideal.
- `util`: auxilary code shared among libraries.

The sources for the modules are in [`src/`](src).  Each module is
structured as follows:

```
SQIsign
└── src
    └── <module_name>
        ├── <arch>
        │   ├── generic
        │   └── lvl1
        ├── opt
        │   ├── generic
        │   ├── lvl1
        │   ├── lvl1_var1
        │   ├── lvl3
        │   └── lvl5
        └── ref
            └── generic
```

where:

- `<module_name>` is the name of the module.
- `<arch>` are optional architecture-specific implementations of the
  module (e.g., `broadwell` for code using assembly instructions
  specific to the Broadwell platform).
- `opt` and `ref` are the portable *optimized* and *reference*
  implementations.
- `lvl1`, `lvl3`, `lvl5` are parameter-dependent implementations of
  the module, corresponding to NIST levels 1, 3 and 5, respectively.
- `lvl1_var1` is a variant of `lvl1`, e.g., using a different prime
  characteristic. The naming is free, and implementors are encouraged
  to choose more explicit naming, e.g., `lvl1_varp6983` for the
  variant using the `p6983` prime defined in the SQIsign AC20 variant.
- `generic` is a parameter-independent implementation of the module.
  If no folder is named like the currently selected variant (see
  [Build](README.md#Build)), then this is compiled instead.
  
Each of the folders above is allowed to be a symlink. E.g., if a
module has no separate optimized and reference implementation, then
`opt` can be a symlink to `ref`. Other example: a module's code only
depends on the field size, but not the specific prime, then
`lvl1_varp6983` could be a symlink to `lvl1`.

### Contents of a module

The leaf folders described above should arrange code as described
below.  We use the `generic` implementation of the `uintbig` module as
an example.

```
generic
├── bench
│   ├── CmakeLists.txt
│   ├── bench1.c
│   └── bench2.c
├── include
│   └── uintbig.h
├── test 
│   ├── CmakeLists.txt
│   ├── test1.c
│   └── test2.c
├── CmakeLists.txt
├── internal_header.h
├── soruce1.c
└── soruce2.c
```

where:

- `include/` shall contain a **unique header file** named
  `<module_name>.h`, where `<module_name>` is the name of the module.
  This header contains the public API of the module, and is the only
  header that can be included by other modules (e.g., via `#include
  <uintbig.h>`). These files must contain extensive doxygen-formatted
  documentation describing the module, see
  [Documentation](#Documentation).
- `bench` and `test` contain one executable per file, containing,
  well, benchmarks and unit tests. Refer to [Benchmarks](#Benchmarks)
  and [Tests](#Tests) for instructions on how to write these.
- Internal headers for the private use of the module, such as
  `internal_header.h` go to the root. Include these using `#include
  "internal_header.h"`.
- The implementation of the module also goes into the root.


## Tests

It is important to have extensive test coverage of the whole software.
Each module must have its own unit tests, as well as integration tests
to ensure consistency across the modules.

### Unit tests

These go into `src/<module_name>/<ref|opt|...>/<generic|lvl1|...>/test/`.
Refer to ... for an example of how to write tests.

### Integration tests

These go into `test/`.  Refer to
[`test/test_sqisign.c`](test/test_sqisign.c) for an example.

### Known Answer Tests (KAT)

KATs help validate consistency across implementations. By ensuring
that, e.g., the optimized and reference implementation produce the
same signatures.

See [Known Answer Tests in README.md](README.md#Known Answer Tests (KAT)).

## Benchmarks

Benchmarks for a module go into
`src/<module_name>/<ref|opt|...>/<generic|lvl1|...>/bench/`.  Global
benchmarks go...

## Documentation

Use [Doxygen headers](https://www.doxygen.nl/manual/docblocks.html)
for documentation.

All code should be extensively documented.  The public module headers
**MUST** be thoroughly documented.

CI automatically builds a PDF of the doc every time code is pushed.
To download the PDF, go to
[Actions](https://github.com/SQIsign/sqisign-nist/actions), click on
the workflow run you're interested in, then go to Artifacts -> docs
(see figure).

![](https://user-images.githubusercontent.com/149199/231756751-0f2780f8-33fe-4db9-8800-b5f145423b65.png)

## Branches and pull requests

Always work on topic branches, never push work in progress on the
`main` branch.  Once a task / issue / work unit is completed, create a
pull-request and ask your team leader for a review.

## Coding style

- **C version**: All code must compile cleanly as *C99*, without
  emitting any warnings, using recent versions of GCC and clang.

- **Names**: Externally visible functions and types should be prefixed
  with the name of the module they belong to.

- **Aliases**: Do use `typedef` with descriptive names whenever it
  makes any nonzero amount of sense to do so.
  Avoid `typedef`s for things other than `struct`s (or elementary data
  types); they have a tendency to break for array and pointer types if
  programmers are not aware of the nature of the underlying type.

- **Parameters**: Output arguments, if any, should always come first.
  Input arguments should generally be marked `const`. Objects of types
  which typically fit into registers should be passed and returned by
  value, larger objects by reference (i.e., as a pointer).
  If certain arguments often appear together, it may be an indication
  that they should be wrapped as a `struct`.

- **Global variables**: Global *constants* are acceptable if needed,
  especially within modules whose code already implicitly relies on
  the same constants anyway (primary example: 𝔽ₚ). It is often a good
  idea to group global constants in a meaningful `struct` and write
  the code such that the struct could easily be replaced by a runtime
  variable at a later point.
  Global *state* (modifiable global variables), on the other hand, is
  strictly forbidden.

- **Whitespace**: Try not to mix tabs and spaces. Line endings
  should be UNIX-style (i.e., `\n` rather than `\r\n`). Whitespace
  characters at the end of a line, or by themselves on an otherwise
  empty line, are to be avoided.
