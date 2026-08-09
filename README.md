# sparrowhawk-amr
Offline k-mer detector of antimicrobial resistance genes for the sparrowhawk toolkit, based on the [AMRFinderPlus](https://github.com/ncbi/amr) database and written in Rust.

---
## Disclaimer :warning: :construction:
This is a **work in progress** project. This in particular implies:

- Not all the main features we want are yet implemented.
- Code might be messy, and not even documented.
- General documentation on how to install and use the tool might be short or even missing.
- Finding unexpected errors/behaviour or bugs should not be a surprise.

These (and potentially other) items will be progressively fixed before version 1.0.

---

## sparrowhawk?
Sparrowhawk was at one time the Archmage of [Earthsea](https://en.wikipedia.org/wiki/Earthsea).
Also, the [sparrowhawk](https://en.wikipedia.org/wiki/Eurasian_sparrowhawk) (*Accipiter nisus*) is a bird of prey native to Europe (and the island of Gont).

# Description

**Note:** this repository contains the antimicrobial resistance (AMR) detector of the sparrowhawk toolkit. If you are looking for the Rust-based genomic assembler, see [sparrowhawk-asm](https://github.com/bacpop/sparrowhawk-asm); for the web implementation of the toolkit, see [sparrowhawk](https://github.com/bacpop/sparrowhawk).

sparrowhawk-amr detects AMR (as well as stress and virulence) genes in assembled bacterial genomes by k-mer matching against the [NCBI AMRFinderPlus database](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/), completely offline: once the database is fetched and indexed, no call to AMRFinderPlus (or BLAST/HMMER) is needed, so the detector can run natively or in a browser.

Current **main features**:
- Fetching of the AMRFinderPlus reference database from the NCBI FTP (`db fetch`).
- Construction of a compact, binary k-mer index from it (`index build`), with a DNA or protein alphabet (k=31 by default), and a selectable subset of reference types (AMR, stress, and/or virulence), as well as index inspection commands (`index stats`, `index report-map`, `index unit-stats`).
- Detection over the contigs directly (`detect direct`), or over genes called from them (`detect cds`) with [orphos](https://github.com/vrbouza/orphos) (which can also be run standalone with `genes call`), including protein-level detection. Results are reported as JSON. Note that no point-mutations are supported currently.
- Evaluation and debugging tooling against native AMRFinderPlus runs (`eval` subcommands), plus a Python-based benchmark suite in [`benchmark`](./benchmark) (see its own README).
- Compilation both to native and WebAssembly targets: the wasm build exposes an `AmrDetector` interface to load an index and detect from the browser.

# Installation
Currently the only option is to compile from source. You will need the [rust toolchain](https://www.rust-lang.org/tools/install) installed in your system. Development has been done only on x86_64 GNU/Linux-based systems, and most surely will probably stay that way (i.e. no other systems have been tested). To get an in principle working version (up to some degree), always clone a versioned tag. E.g. to get version vX.Y.Z, you could use:

```
git clone --branch vX.Y.Z https://github.com/bacpop/sparrowhawk-amr.git
cd sparrowhawk-amr
cargo build --release
```

This should place your compiled binary inside `target/release`. As with the rest of the toolkit, you can also compile it to the WebAssembly target `wasm32-unknown-unknown` (the wasm-specific code paths are selected automatically through `cfg(target_family = "wasm")`):

```
rustup target add wasm32-unknown-unknown
cargo build --release --target wasm32-unknown-unknown
```

# Usage
You can see the available subcommands and their arguments with

```
./sparrowhawk-amr --help
```

(as well as with `--help` on each subcommand). A minimal end-to-end run — fetch the database, index it, and detect on an assembly — could be:

```
./sparrowhawk-amr db fetch --out-dir ./amrfinder_db
./sparrowhawk-amr index build --db-dir ./amrfinder_db --out ./amr_index.bin
./sparrowhawk-amr detect direct --index ./amr_index.bin --fasta ./contigs.fasta --sample-name mysample
```

The detections are printed to stdout as JSON. To detect on called genes instead of raw contigs, use `detect cds` with the `--assembly` argument (gene calling is done internally with orphos), and to compare against a native AMRFinderPlus run, check the `eval` subcommands and the [benchmark suite](./benchmark).
