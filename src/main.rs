#![allow(warnings)]
use cherry_red::network::Network;
use std::fs;
use std::env;
fn main() {
    let args: Vec<String> = env::args().collect();
    //     ---------        set up environment variables
    let mut is_debug = false;
    let mut is_exact = false;
    let is_gen_mode;
    let mut is_rank = false;
    let mut rank_t_opt: Option<f32> = Some(0.5);
    let mut req_args: Vec<&str> = vec![];
    //       --------------             parse command
    match args[1].as_str() {
        "gen" => is_gen_mode = true,
        "read" => is_gen_mode = false, 
        _ => {help("error: provide command for mode.");return},
    }
    //    -------------          parse option flags
    //short form and long form
    for (i,arg) in args.iter().enumerate() {
        if i == 0 || i == 1 {
            continue;
        }
        let mut v: Vec<&str> = arg.split('=').collect();
        v[0] = v[0].trim_start_matches('-');
        if arg.starts_with("--") {
            match v[0] {
                "debug" => is_debug = true,
                "verbose" => is_debug = true,
                "exact" => is_exact = true,
                "rank" => {is_rank = true;
                        if let Some(e)= v.get(1) {
                            rank_t_opt = Some(e.parse::<f32>().unwrap());
                        }},
                "help" => {help("printing help...");return},
                _ => {help("unknown flag encountered"); return},
            }
        } else if arg.starts_with('-') {
            for c in v[0].chars() {
                match c {
                    'd' => is_debug = true,
                    'v' => is_debug = true,//legacy "verbose"
                    'e' => is_exact = true,
                    'r' => {is_rank = true;
                        if let Some(e)= v.get(1) {
                            rank_t_opt = Some(e.parse::<f32>().unwrap());
                        }},
                    'h' => {help("printing help...");return},
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
        if is_rank {
            call_gen_mode(is_debug, is_exact, rank_t_opt, leaves, reticulations, distance);
        } else {
            call_gen_mode(is_debug, is_exact, None, leaves, reticulations, distance);
        }
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
        // ---------   read    call distance with selected options and arguments
        if is_rank {
            call_read_mode(is_debug, rank_t_opt, newick1, newick2);
        } else {
            call_read_mode(is_debug, None, newick1, newick2);
        }
    }
}
fn help(message: &str) {
    println!("exit message: {}",message);
    println!("
    Cherry distance program calculates and outputs cherry distance to stdout. It runs in two modes: 

    COMMANDS
    gen     random generation mode
    read     extended format newick string parsing mode
    
    Random generation mode will generate a network with the specified number of leaves and up to the specified number of reticulations (see --exact option below), then randomly modifies the generated network by the given distance. The calculated distance may be less than given (see explanation in manual), though only observed in very small network sizes. 

    Extended format newick string parsing mode takes two required arguments, two files which should each contain a modified extended Newick format string. The modified newick format should not include edge lengths or hybrid nodes types.  

    The program calculates cherry distance exactly by default, however the ranking heuristic can be activated to see a speedup in calculation at the cost of some accuracy. The threshold level is set at 0.5 by default.

    USAGE
    1. random reneration mode
    usage: gen [options] <leaves> <reticulations> <distance> 

    2. extended newick format string parsing mode
    usage: read [options] <file1> <file2>
    
    OPTIONS
    -r     Rank, any mode, which specifies that the ranking heuristic mode should be used. The program calculates cherry distance exactly by default, this is the optional, inexact, faster calculation. 
    -d     Debug, `gen` mode, prints out the generated random networks.
    -e     Exact, `gen` mode, specifies that the exact number of reticulations requested should be reached in randomly generated network, note this may increase runtime.
    -h     Help, any mode, print this help guide.   

    ");
}
fn call_gen_mode(debug:bool, 
                exact:bool, 
                rank_t_opt: Option<f32>,
                leaves:usize, 
                reticulations:usize, 
                distance:usize) {
    let dist;
    let n1 = Network::new_random(leaves, reticulations,exact);
    let (n1, n2) = Network::random_modify(n1, distance);
    if debug {//only prints when rank is activated
        println!("Network 1: {n1}");
        println!("Network 2: {n2}");
        if let Some(rank_t) = rank_t_opt {
            println!("Ranking threshold value: {}",rank_t);
        }
    }
    dist = Network::find_cherry_distance(n1,n2,rank_t_opt);      
    //last output
    println!("Cherry distance of random networks: {dist}");
}
fn call_read_mode(debug:bool, rank_t_opt:Option<f32>, new1: String, new2: String) {
    let n1 = Network::parse_newick(&new1);
    let n2 = Network::parse_newick(&new2);
    //calculate distance
    if debug {
        println!("Parsed network 1: {n1}");
        println!("Parsed network 2: {n2}");
        if let Some(rank_t) = rank_t_opt {//only prints when rank is activated
            println!("Ranking threshold value: {}",rank_t);
        }
    }
    let dist = Network::find_cherry_distance(n1,n2,rank_t_opt);
    //last output
    println!("Cherry distance of newick networks: {dist}");
}