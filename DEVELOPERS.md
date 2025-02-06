# Developer guidelines

Please read carefully before contributing to this repo.

## Code structure

The source code is in the [`src/`](src) directory and it is split into the
following modules:
- `common`: common code for AES, SHAKE, (P)RNG, memory handling. Every
  module that needs a hash function, seed expansion,
  deterministic alea for tests, should call to this module.
- `mp`: code for saturated-representation multiprecision arithmetic.
- `gf`: GF(p^2) and GF(p) arithmetic.
- `ec`: elliptic curves, isogenies and pairings. Everything that is
  purely finite-fieldy.
- `precomp`: constants and precomputed values.
- `quaternion`: quaternion orders and ideals.
- `hd`: code to compute (2,2)-isogenies in the theta model.
- `id2iso`: code for Iso <-> Ideal.
- `verification`: code for the verification protocol.
- `signature`: code for the key generation and signature protocols.

### Contents of a module

Each module is comprised of *implementation types* and common code. An
implementation type refers to the *portable optimized*, *portable reference* or
any architecture-specific implementation of the module; a module must contain at
least one implementation type. Each implementation type must be in its own
directory within the directory of the module. The optimized and reference
implementation types must be placed in the `opt` and `ref` directories,
respectively; there is no rule for naming other architecture-specific
implementations. Common code refers to optional generic code that is shared
among all implementation types. Common code is placed in special directories
within the directory of a module: header files in the `include` directory and
source files in the `<module_name>x` directory, where `<module_name>` is the
name of the module. An example of a module is given below:
```
src
└── <module_name>
    ├── include
    ├── <module_name>x
    ├── opt
    ├── ref
    └── <arch>
```
where:
- `<module_name>` is the name of the module.
- `opt` and `ref` are the portable optimized and reference
  implementation types, respectively.
- `<arch>` is an optional architecture-specific implementation type of the
  module (e.g., `broadwell` for code using assembly instructions
  specific to the Broadwell platform).
- `include` contains header files common to all implementation types.
- `<module_name>x` contains source files common to all implementation types
  (i.e., `opt`, `ref` and `<arch>`).

Header files in the `include` directory above can be included by other modules
and must contain extensive doxygen-formatted documentation describing the
functions declared there; see [Documentation](#Documentation). Any
implementation-type directory above is allowed to be a symlink; e.g., if a
module has no separate optimized and reference implementation, then
`opt` can be a symlink to `ref`.

Similar to a module, each implementation type is comprised of implementation 
*variants* and common code. A variant refers to either a *generic*
implementation, an implementation whose parameters are defined by one of the
NIST levels (i.e., 1, 3 or 5) or a variation of the latter. An implementation
type must contain at least one variant. Each variant must be in its own
directory within that of the implementation type. The generic variant must be
placed in the `generic` directory and variants corresponding to NIST levels 1,
3 and 5 are placed in the directories `lvl1`, `lvl3` and `lvl5`, respectively;
there is no rule for naming the directory of a NIST variation, but
implementors are encouraged to choose informative namings. Common code refers to
optional variant-independent code that is shared among all variants of the same
implementation type. Common code is placed in special directories within that of
the implementation type: header files in the `include` directory and source
files in the `lvlx` directory. Expanding on the example above, we show the
details of its implementation types:
```
src
└── <module_name>
    ├── include
    ├── <module_name>x
    ├── opt
    │   ├── include
    │   ├── lvlx
    │   ├── lvl1
    │   ├── lvl1_var1
    │   ├── lvl3
    │   └── lvl5
    ├── ref
    │   ├── generic
    │   └── lvl1
    └── <arch>
        ├── include
        ├── lvlx
        ├── lvl3
        └── lvl5
```
where:
- `lvl1`, `lvl3`, `lvl5` are implementations of NIST levels 1, 3 and 5,
  respectively, for the corresponding implementation type.
- `lvl1_var1` is a variation of `lvl1` for the `opt` implementation type (e.g.,
  using a different prime characteristic).
- `opt/include` contains header files common to all variants in the `opt`
  implementation type (i.e., `lvl1`, `lvl1_var1`, `lvl3` and `lvl5`).
  Similarly, `<arch>/include` for all variants in the `<arch>` implementation
  type (i.e., `lvl3` and `lvl5`).
- `opt/lvlx` contains source files common to all variants in the `opt`
  implementation type. Similarly, `<arch>/lvlx` for all variants in the `<arch>`
  implementation type.
- `generic` contains a parameter-independent implementation of the `ref`
  implementation type.

As the name suggests, the `generic` variant is a generic implementation which
does not depend on the parameters defined by the NIST levels or any variation
of these. If this directory is present, all other parameter-dependent
implementations are ignored and the `generic` implementation is built instead.
As with modules, header files in the `include` directory of an implementation
type (e.g., `opt/include` and `<arch>/include` above) can be included by other
modules and must contain extensive doxygen-formatted documentation describing
the functions declared there.

Each implementation variant must be organized as follows:
- Header files that can be included by other modules are placed in the `include`
  directory. These files must contain extensive doxygen-formatted documentation
  describing the functions declared there.
- Source files of the implementation and their private internal header files are
  placed directly in the implementation variant directory.
- Source files of unit tests and their private internal header files are placed
  in the `test` directory. Refer to [Tests](#Tests) for instructions on how to
  write these.

Common code (in `lvlx`) for all variants in an implementation type follows the
same organization as above, with the exception that `lvlx` never contains an
`include` directory. This role is taken by the `include` directory in the
implementation type. Below is an example with the detailed organization of the
common code and the `lvl1` variant for the `ref` implementation type of a
module:
```
<module_name>
├──ref
│  ├── include
│  │   └── header_ref.h
│  ├──lvlx
│  │  ├── test
│  │  │   ├── test_internal_header_ref.h
│  │  │   │   ...
│  │  │   ├── test1_ref.c
│  │  │   └── test2_ref.c
│  │  ├── internal_header_ref.h
│  │  ├── source1_ref.c
│  │  └── source2_ref.c
│  ├──lvl1
│  │  ├── include
│  │  │   └── header_ref_lvl1.h
│  │  ├── test
│  │  │   ├── test_internal_header_ref_lvl1.h
│  │  │   │   ...
│  │  │   ├── test1_ref_lvl1.c
│  │  │   └── test2_ref_lvl1.c
│  │  ├── internal_header_ref_lvl1.h
│  │  ├── source1_ref_lvl1.c
│  │  └── source2_ref_lvl1.c
│  ├──lvl3
│  └──lvl5
```

Finally, common code for a module must be organized as follows:
- Header files that can be included by other modules are placed in the `include`
  directory. As mentionde before, these files must contain extensive
  doxygen-formatted documentation describing the functions declared there.
- Source files and their private internal header files are placed in the 
  `<module_name>x` directory.
- Source files of unit tests and their private internal header files are placed
  in the `<module_name>x/test` directory. Again, refer to [Tests](#Tests) for
  instructions on how to write these.

The example below shows the detailed organization of the common code of a
module:
```
<module_name>
├── include
│   └── header.h
├── <module_name>x
│   ├── test
│   │   ├── test_internal_header.h
│   │   │   ...
│   │   ├── test1.c
│   │   └── test2.c
│   ├── internal_header.h
│   ├── source1.c
│   └── source2.c
├── opt
└── ref
```

## Tests

It is important to have extensive test coverage of the whole software.
Each module must have its own unit tests, as well as integration tests
to ensure consistency across the modules.

### Unit tests

These go in the `src/<module_name>/<module_name>x/test` and 
`src/<module_name>/<ref|opt|...>/<generic|lvlx|lvl1|...>/test/` directories.
Refer to [`src/gf/gfx/test/test_fp.c`](src/gf/gfx/test/test_fp.c) for an example
of how to write tests.

### Integration tests

These go in the `test/` directory. Refer to
[`test/test_sqisign.c`](test/test_sqisign.c) for an example.

### Known Answer Tests (KAT)

KATs help validate consistency across implementations. By ensuring
that, e.g., the optimized and reference implementation produce the
same signatures. KATs are generated by executing `PQCgenKAT_sign_<level>` in
the `apps` directory. KAT tests go in the `test/` directory.

## Benchmarks

Benchmarks for a module go in the same directories as for tests.
Global benchmarks go in the `apps` directory; e.g.,
[`apps/benchmark.c`](apps/benchmark.c).

## Documentation

Use [Doxygen headers](https://www.doxygen.nl/manual/docblocks.html)
for documentation.

All code should be extensively documented.  The public module headers
**MUST** be thoroughly documented.

CI automatically builds a PDF of the doc every time code is pushed.
To download the PDF, go to
[Actions](https://github.com/SQIsign/sqisign-nist2/actions), click on
the workflow run you're interested in, then go to Artifacts -> docs
(see figure). PDFs are retained for 2 days.

![](https://user-images.githubusercontent.com/149199/231756751-0f2780f8-33fe-4db9-8800-b5f145423b65.png)

## Branches and pull requests

Always work on topic branches, never push work in progress on the
`main` branch.  Once a task / issue / work unit is completed, create a
pull-request and ask for at least one review.

## Coding style

- **C version**: All code must compile cleanly as *C11*, without
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

- **Formatting**: This project uses
  [`clang-format`](https://clang.llvm.org/docs/ClangFormat.html) to
  format the code.  From the root of the project run the following
  command:
  
  ```
  find ./src -path ./src/precomp -prune -type f -o -iname '*.h' -o -iname '*.c' | xargs clang-format -i 
  ```
  
  to automatically format all appropriate files with `clang-format`.
  
  If you want, you can install a [pre-commit
  hook](https://pre-commit.com/) to ensure that your work is correctly
  formatted before pushing
  
  ```
  pre-commit install
  ```
  
  Will use the `.pre-commit-config.yaml` file.
