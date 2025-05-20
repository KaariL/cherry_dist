# cherry-red

**Cherry Red**(uced) is a program to calculate cherry distance metric between two binary, level-1, orchard networks. The paper detailing the implementation and use of this code is:

> *Fast calculation of cherry distance on level-1 orchard networks: optimization, heuristic and implementation*. Kaari Landry and Olivier Tremblay-Savard. 2025.
> 
## Installation

### With the precompiled binary directly
```
git clone https://github.com/KaariL/cherry_dist
/.cherry-red <usage>
```
### With Cargo
[Install Cargo](https://doc.rust-lang.org/cargo/getting-started/installation.html) and use either `run`, or `build`:
- `cargo run <usage>` or,
- `cargo build` (use options like `--release` to optimize build).
#### Manual
With `cargo` installed, run `cargo docs --open` to see documentation for modules and public functions.
## Usage
RedCherry calculates and outputs cherry distance to stdout. It runs in two modes, details below, with the following usage:

1.  `gen` random generation mode 
	`gen [options] <leaves> <reticulations> <distance>`
2. `read` extended format Newick string parsing mode
	`read [options] <file1> <file2>`

The available options are:
`-r,--rank` Rank, any mode, which specifies that the ranking heuristic mode should be used. The program calculates cherry distance exactly by default, this is the optional, inexact, faster calculation. Choose a threshold in (0,1), set to 0.5 by default. Examples: `-r=0.75`, `--rank=0.25`
`-d, --debug` Debug, any mode, prints out the generated or parsed networks.
`-e, --exact` Exact, `gen` mode, specifies that the exact number of reticulations requested should be reached in randomly generated network, note this may increase runtime.
`-h, --help` Help, any mode, print this help guide.

### Extended format Newick string parsing mode
[Extended format Newick](http://www.biomedcentral.com/1471-2105/9/532) string parsing mode takes two required arguments, two files which should each contain a modified extended Newick format string. The modified Newick format should not include edge lengths or hybrid nodes types.
### Random generation mode
Random generation mode will generate a network with the specified number of leaves and up to the specified number of reticulations (see --exact option below), then randomly modifies the generated network by the given distance. The calculated distance may be less than given (see explanation in manual), though only observed in very small network sizes.

## Testing example  
There are 3 testing files available in `\examples` transcribed from networks found in [this work](https://doi.org/10.1111/jipb.13246).
Running the compiled binary directly:
```
./cherry-red read -r=0.1 examples/n1.txt examples/n2.txt
```

Using `cargo run`:
```
cargo run read examples/n0.txt examples/n0.txt
# Cherry distance of newick networks: 0
```

Using `cargo build`:
```
cargo build --release
cd target/release
./cherry-red read ../../examples/n1.txt ../../examples/n2.txt
# Cherry distance of newick networks: 39
```

---
