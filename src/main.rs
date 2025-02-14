use cherry_dist::network::Network;
use std::fs;
use std::env;

fn main() {
    env::set_var("RUST_BACKTRACE", "1");
    let args: Vec<String> = env::args().collect();
    //     ---------        set up environment variables
    let mut is_debug = false;
    let mut is_exact = false;
    let mut is_gen_mode = false;
    let mut req_args: Vec<&str> = vec![];
    //       --------------             parse command
    match args[1].as_str() {
        "gen" => is_gen_mode = true,
        "new" => is_gen_mode = false, 
        _ => {help("error: provide command for mode.");return},
    }
    //    -------------          parse option flags
    //short form and long form, in alphabetical order
    for (i,arg) in args.iter().enumerate() {
        if i == 0 || i == 1 {
            continue;
        }
        if arg.starts_with("--") {
            let option = arg.as_str().trim_start_matches('-');
            println!("option: {option}");
        } else if arg.starts_with('-') {
            let options = arg.trim_start_matches('-');
            for c in options.chars() {
                match c {
                    'd' => is_debug = true,
                    'v' => is_debug = true,//legacy "verbose"
                    'e' => is_exact = true,
                    'h' => {help("printing help...");return},
                    'm' => {println!("a link to the manual.");},//TODO
                    _ => {print!("{c} ");help("unknown flag encountered"); return},
                }
            }
        } else {
            req_args.push(arg);
        }
    }
    //      -----------           parse arguments and call
    if is_gen_mode {
        let mut leaves = 0;
        let mut reticulations = 0;
        let mut distance = 0;
        for (i,req_arg_result) in req_args.into_iter().enumerate() {
            if let Ok(req_arg) = req_arg_result.parse::<usize>() {
                match i {
                    0 => leaves = req_arg,
                    1 => reticulations = req_arg,
                    2 => distance = req_arg,
                    _ => {help("Not the required number of arguments for random generation mode"); return},
                }
            } else {
                help("gen mode required argument is not numeric."); return
            }
        }
        //      ------------   gen   call distance with given args and options 
        call_gen_mode(is_debug, is_exact, leaves, reticulations, distance);
    } else {
        if req_args.len() != 2 {
            help("Not the required number of arguments for extended newick format string parser mode"); return
        }
        let file_result_1 = fs::read_to_string(req_args[0]);
        let newick1 = match file_result_1 {
            Ok(new_string) => new_string,
            Err(_) => {help("Could not read file 1."); return}
        };
        let file_result_2 = fs::read_to_string(req_args[1]);
        let newick2 = match file_result_2 {
            Ok(new_string) => new_string,
            Err(_) => {help("Could not read file 2."); return}
        };
        // ---------   new     call distance with selected options and arguments
        call_new_mode(is_debug, newick1, newick2);
    }
}
fn help(message: &str) {
    println!("exit message: {}",message);
    println!("
    Cherry distance program calculates and outputs cherry distance to stdout. It runs in two modes: 

    COMMANDS
    gen     random generation mode
    new     extended format newick string parsing mode
    
    Random generation mode will generate a network with the specified number of leaves and up to the specified number of reticulations (see --exact option below), then randomly modifies the generated network by the given distance. The calculated distance may be exactly 2 less than given (see manual elaboration). 

    Extended format newick string parsing mode takes two required arguments, two files which should each contain an extended Newick format string. The program will parse each into a network, outputting the resulting parsed network each into an output file. 

    USAGE
    1. random reneration mode
    usage: ./cherry-dist gen [options] <leaves> <reticulations> <distance> 

    2. extended newick format string parsing mode
    usage: ./cherry-dist new [options] <file1> <file2>
    
    OPTIONS
    -d, --debug     verbose printing mode outputs mainly runtime information
    -e, --exact     option only available on random generation mode, specifies that the exact number of reticulations requested should be reached in randomly generated network, note this may increase runtime
    -h, --help      print this help guide   
    -m, --manual    links to version of manual on web

    ");
}
fn call_gen_mode(debug:bool, 
                exact:bool, 
                leaves:usize, 
                reticulations:usize, 
                distance:usize) {
    let dist;
    let n1 = Network::new_random(leaves, reticulations,exact);
    let (n1, n2) = Network::random_modify(n1, distance);
    if debug {
        println!("Network 1: {n1}");
        println!("Network 2: {n2}");
    }
    dist = Network::find_cherry_distance(n1,n2);      
    //last output1
    println!("Cherry distance of random networks: {dist}");
}
fn call_new_mode(debug:bool, new1: String, new2: String) {
    let n1 = Network::parse_newick(&new1);
    let n2 = Network::parse_newick(&new2);
    //calculate distance
    let dist = Network::find_cherry_distance(n1,n2);
    //last output
    println!("Cherry distance of newick networks: {dist}");
}