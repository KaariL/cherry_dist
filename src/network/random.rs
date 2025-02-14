//! Random Submodule (of Network )
//!
//! Generates a random network topology that is level-1 and orchard
//! 
//! `random` submodule contains all functionality associated with generating 
//! random trees and for attaching additional edges to form reticulations.
//! 

use rand::{distributions::Uniform, Rng};
use rand::thread_rng;//for shuffle
use rand::seq::SliceRandom;//for shuffle
use std::collections::BTreeSet;
use std::cmp;
//
use crate::network::Network;
use crate::network::Node;
use crate::network::gen_label;
use crate::network::distance::DTree;

// --------------------------- PUBLIC
pub fn new_random(size: usize, r: usize, exact: bool) -> Network<Node<String>> {
    const MAX_TRY_FOR_EXACT: i32 = 100;
    if size < 2 {panic!("cannot construct network with less than 2 leaves")}
    //NOTE:check r isnt too big
    //size is # leaves
    //v is # nodes
    let v = 2*size-1;
    //generate the random numbers needed
    let mut rands: Vec<usize> = vec![]; 
    let mut rng = rand::thread_rng();
    let dist = Uniform::from(1usize..100usize);//[1,99] (percent)
    for _ in 0..v {
        rands.push(rng.sample(dist));
    } 
    //------------------------build tree and add reticulations
    let mut result: Network<Node<String>> = Network::new(size);
    if exact {
        //build tree
        for i in 0..MAX_TRY_FOR_EXACT {
            if i == MAX_TRY_FOR_EXACT-1 {
                panic!("could not add exact number of reticulations in {} tries", MAX_TRY_FOR_EXACT);
            }
            let mut net = Network::new(size);
            net = new_random_tree(&rands, net, 0, size, 0);
            net.update_tree_clusters();
            //add rets
            result = add_random_rets(net, r);
            if result.get_ret_nodes().len() == r {
                break
            }
        }
    } else {
        //build tree
        let mut net = new_random_tree(&rands, result, 0, size, 0);
        net.update_tree_clusters();
        //add rets
        result = add_random_rets(net, r);
    }
    result
}
pub fn new_random_tree_from_x(size: usize, mut x: Vec<String>) -> Network<Node<String>> {
    if size < 2 || x.len() < 2 {
        panic!("cannot construct network with less than 2 leaves")
    }
    if size != x.len() {
        panic!("num leaves and num labels do not match")
    }
    let v = 2*size-1;
    //generate the random numbers needed
    let mut rands: Vec<usize> = vec![]; 
    let mut rng = rand::thread_rng();
    let dist = Uniform::from(1usize..100usize);//[1,99] (percent)
    //TODO this is a bad way to do it. not granular enough
    for _ in 0..v {
        rands.push(rng.sample(dist));
    } 
    //------------------------build tree
    let mut result: Network<Node<String>> = Network::new(size);
    //build tree
    x.shuffle(&mut thread_rng());
    let mut tree = new_random_tree_from_x_helper(&rands, &x, result, 0, size, 0);
    tree.update_tree_clusters();
    tree
}
pub fn new_random_from_level_seq(seq: Vec<u16>, l: usize) -> Network<Node<String>> {
    //error check: 
    if l < 4 * (seq.len()-1) + 2 {panic!("not enough leaves");}
    let dtree = DTree::build_from_level_seq(seq);
    let mut base:Network<Node<String>> = base_from_dtree(dtree);
    grow(&mut base,l);
    base.update_tree_clusters();
    base
}
// --------------------------- building the tree
fn new_random_tree_from_x_helper(
    cur_rands: &[usize],
    cur_x: &[String],
    mut cur_net: Network<Node<String>>,
    cur_index: usize,
    cur_size: usize,//# of leaves
    par_i: usize,
) -> Network<Node<String>> {
    //assumes x is already shuffled
    let mut rng = rand::thread_rng();
    let dist = Uniform::from(0usize..usize::MAX);
    //assumes Network will not be singleton
    if cur_size < 2 {//BASE CASE
        //create leaf and fill cur_node position in Network
        //choose from x and remove it
        let rand_i = rng.sample(dist);
        let rand_label = cur_x[0].clone();
        let ( _ , cur_x) = cur_x.split_at(1);//remove label once used
        cur_net.add_leaf(cur_index, rand_label);
        cur_net.add_edge(par_i, cur_index);
    } else {
        //create+add internal node at cur_node position in Network
        if cur_index == 0 {
            cur_net.add_node_i(0);
        } else {
            cur_net.add_node_i(cur_index);
            cur_net.add_edge(par_i, cur_index);
        }
        //   ---------   recursive calls for left and right subnetworks
        //determine size (# leaves) of subNetworks
        let left_pc = cur_rands[0];
        let left_size = cmp::max(1, left_pc*cur_size/100);
        let right_size = cur_size - left_size;
        let ( _ , cur_rands) = cur_rands.split_at(1);
        //left and right of rands
        let (left_rands, right_rands) = cur_rands.split_at(2*left_size-1);
        //left and right indices
        let left_index = cur_index + 1;
        let right_index = cur_index + (2*left_size-1) + 1;
        //left and right of x 
        let (left_x, right_x)= cur_x.split_at(left_size);
        //recursive calls
        cur_net = new_random_tree_from_x_helper(left_rands, left_x, cur_net, left_index, left_size, cur_index);
        cur_net = new_random_tree_from_x_helper(right_rands, right_x, cur_net, right_index, right_size, cur_index);
    }
    cur_net
}
fn new_random_tree(
        cur_rands: &[usize],
        mut cur_net: Network<Node<String>>,
        cur_index: usize,
        cur_size: usize,//# of leaves
        par_i: usize,
    ) -> Network<Node<String>> {
    //assumes Network will not be singleton
    if cur_size < 2 {//BASE CASE
        //create leaf and fill cur_node position in Network
        cur_net.add_leaf(cur_index, gen_label());
        cur_net.add_edge(par_i, cur_index);
    } else {
        //create+add internal node at cur_node position in Network
        if cur_index == 0 {
            cur_net.add_node_i(0);
        } else {
            cur_net.add_node_i(cur_index);
            cur_net.add_edge(par_i, cur_index);
        }
        //   ---------   recursive calls for left and right subnetworks
        //determine size (# leaves) of subNetworks
        let left_pc = cur_rands[0];
        let left_size = cmp::max(1, left_pc*cur_size/100);
        let right_size = cur_size - left_size;
        let ( _ , cur_rands) = cur_rands.split_at(1);
        //left and right of rands
        let (left_rands, right_rands) = cur_rands.split_at(2*left_size-1);
        //left and right indices
        let left_index = cur_index + 1;
        let right_index = cur_index + (2*left_size-1) + 1;
        //recursive calls
        cur_net = new_random_tree(left_rands, cur_net, left_index, left_size, cur_index);
        cur_net = new_random_tree(right_rands, cur_net, right_index, right_size, cur_index);
    }
    cur_net
}
fn add_random_rets(net: Network<Node<String>>, r:usize) -> Network<Node<String>> {
    //only guarunteed to return up to r edge pairs, not exactly r
    if !net.is_tree() { panic!("reticulations cannot be added piecemeal."); }
    let potential_bicons:Vec<Vec<(usize,usize)>> = vec![net.get_edges()];
    if potential_bicons[0].len() < 2 { net } else {
        add_random_rets_helper(net,r, potential_bicons)
    }
}
fn add_random_rets_helper(
                    mut net: Network<Node<String>>, 
                    r: usize, 
                    mut potential_bicons: Vec<Vec<(usize,usize)>>, 
    ) -> Network<Node<String>> {
    if r==0 || potential_bicons.is_empty()  {
        return net;
    }
    //set up random distribution
    let mut rng = rand::thread_rng();
    let dist = Uniform::from(0usize..usize::MAX);
    //b = total across all groups. choose in proportion to group size
    let a = rng.sample(dist);
    let b:usize = potential_bicons.iter().map(|l| l.len()).sum();
    //set up edge selection
    let e1:(usize,usize);let e2:(usize,usize);
    let mut cur_group; let mut g = 0; //g is index of cur group
    let mut edge_meta_i = a % b;
    //start selection
    'edge_1: loop {
        for group in &mut potential_bicons {
            for (i, _) in group.iter().enumerate() {
                if edge_meta_i == 0 {   
                    e1 = group.swap_remove(i);
                    break 'edge_1;
                }
                edge_meta_i -= 1;
            }
            g+=1;
        }
    }
    cur_group = potential_bicons.swap_remove(g);
    //choose e2
    let a = rng.sample(dist);
    let b = cur_group.len();
    let edge_i = a % b;
    e2 = cur_group.swap_remove( edge_i );
    //choose direction of ret edge
    let (_, in_1) = e1; let (_, in_2) = e2; 
    let mut out_edge = e1; let mut in_edge = e2;
    if net.is_below(in_1, in_2) {
        out_edge = e2; in_edge = e1;
    }
    //add the reticulation
    let ret_node_i = net.add_reticulation(out_edge,in_edge,r); 
    let (_, bc) = in_edge;
    let (_, ac) = out_edge;
    cur_group.push((ret_node_i,bc));
    cur_group.push((ret_node_i-1,ac));
    //update potential bicons
    let bicon = net.get_bicon(ret_node_i);
    cur_group.retain(|&(a,b)| !(bicon.contains(&a) && bicon.contains(&b)));
    //var set up
    let mut new_groups:Vec<Vec<(usize,usize)>> = vec![vec![];bicon.len()];
    let mut outer_roots:Vec<Option<usize>> = vec![None;bicon.len()];
    let mut moved = vec![false;cur_group.len()];
    //find "outer roots" group
    for (b_i,v) in bicon.iter().enumerate() {
        //so inefficient
        for (cg_i,(a,b)) in cur_group.iter().enumerate() {
            if a==v {
                outer_roots[b_i] = Some(*b);
                moved[cg_i] = true;//flag for later removal from cur_group
            } 
        }
    }
    for (cg_i, (a,b)) in cur_group.iter().enumerate() {
        //add edges to a new group
        for (r_i,root) in outer_roots.iter().enumerate() {
            if root.is_some() {
                if net.is_below(*a,root.unwrap()) {
                    new_groups[r_i].push((*a,*b));
                    moved[cg_i] = true;
                } else if *b == root.unwrap() {
                    //no reticulation is an outer root so
                    //this is the rooting edge
                    new_groups[r_i].push((*a,*b));
                    moved[cg_i] = true;
                }
            }
        }
    }
    //all other edges are thier own group
    let mut moved_iter = moved.iter();
    cur_group.retain(|_| !(moved_iter.next().unwrap()) );
    new_groups.push(cur_group);
    //add new groups to potential bicons if large enough
    for new_group in new_groups {
        if new_group.len() > 1 {
            //conjecture: 2 connected edges can have ret edge between them
            potential_bicons.push(new_group);
        }
    }
    //update network, down here due to borrowing rules
    net.add_nontriv_bicon(bicon,r);
    //call again for next reticulation
    add_random_rets_helper(net, r-1, potential_bicons)
}
// -------------------- Generate Network from DTree
fn base_from_dtree(dtree: DTree) -> Network<Node<String>> {
    let l = dtree.get_n() * 2;//this is less then n but gets us closer 
    //add root and two leaf children
    let mut base = Network::new(l);
    base.add_node_i(0);
    base.add_anon_leaf_child_of(0,gen_label());
    base.add_anon_leaf_child_of(0,gen_label());
    //processing associates dtree v with l where component is component is added above
    let mut processing:Vec<(usize,usize)> = vec![];
    for c_i in &dtree.get_children_i(0) {
        //just always goes on the left child of Nroot
        let p_1 = base.get_parents_i(1)[0];
        let c_id = bisect_pendant_leaf(&mut base, (p_1,1));
        processing.push((*c_i,c_id));
    }
    while !processing.is_empty() {
        //id in network, i in tree
        let (v_i,l_id) = processing.pop().unwrap();
        //replace v_id with base component
        let r_id = bisect_base_comp(&mut base, l_id);
        for c_i in dtree.get_children_i(v_i) {
            let r_c_id = base.get_children_i(r_id)[0];
            let c_id = bisect_pendant_leaf(&mut base, (r_id,r_c_id));
            processing.push( (c_i,c_id) );
        }
    }
    base
}
//adds a randomly labelled pendant leaf bisecting edge e, returns leaf id
fn bisect_pendant_leaf(n: &mut Network<Node<String>>, e:(usize,usize)) -> usize  {
    //subdivide edge
    //make and add new node p_id
    let p_id = n.add_anon_node();
    n.subdivide_edge(e,p_id);
    //add leaf
    let l_id = n.add_anon_leaf_child_of(p_id,gen_label()); 
    //return the index of the new leaf
    l_id
}
fn bisect_anon(n:&mut Network<Node<String>>, e:(usize,usize)) -> usize{
    //bisect edge with anon internal node, return internal node
    let id = n.add_anon_node();
    n.subdivide_edge(e,id);
    id
}
fn bisect_base_comp(n:&mut Network<Node<String>>, l: usize) -> usize {
    //returns the edge leaving r
    let p_l = n.get_parents_i(l)[0];
    let p_id = bisect_anon(n, (p_l,l));
    //add base root
    let a = n.add_anon_child_of(p_id);
    //add its 2 children
    let b = n.add_anon_child_of(a);
    let c = n.add_anon_child_of(a);
    //make r 
    let r = n.add_anon_child_of(b);
    //connect to r
    n.add_edge(c,r);
    n.add_addn_edge((c,r));
    //add the notriv bicon to the network
    let bicon = vec![a,b,c,r];//not the whole bicon!
    n.add_nontriv_bicon(bicon,r);
    //add 2 anon leaves
    let _ = n.add_anon_leaf_child_of(r,gen_label());
    let _ = n.add_anon_leaf_child_of(b, gen_label());
    let _ = n.add_anon_leaf_child_of(c, gen_label());
    //
    r
}
fn grow(base: &mut Network<Node<String>>, leaves: usize) {
    //calculate how many leaves are already present?
    let n = leaves - 2 - 4*(base.get_ret_nodes().len());
    //set up
    let mut es = base.get_edges();
    let mut rng = rand::thread_rng();
    let dist = Uniform::from(0usize..usize::MAX);//[0,max)
    //processing
    for _ in 0..n {
        //choose random edge
        let a = rng.sample(dist);
        let b = es.len();
        let i = a % b;
        let  e = es.remove(i);
        //add new random leaf to network
        let leaf = bisect_pendant_leaf(base, e);
        //update edges
        let p_l = base.get_parents_i(leaf)[0];
        es.push((p_l,leaf));
        es.push((e.0,p_l));
        es.push((p_l,e.1));
    }
}
