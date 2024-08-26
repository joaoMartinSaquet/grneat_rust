// use std::{collections::TryReserveError, path};


use std::panic::RefUnwindSafe;

use rust_grn::grn::{grn_model, grn_protein};
use rust_grn::evolver::grn_gene::GrnGene;
use plotly::common::{Mode, Title};
use plotly::layout::{Axis, Layout};
use plotly::{traces, Plot, Scatter};
use rand::*;

fn test_grn()
{
    let path: &str = "/home/jmartinsaquet/Documents/code/GRN_rust_java/GRNEAT/SinusExperience/run_24561318703612/grn_49_-27.37311533076093.grn";
    // let test = grn_model::GrnModel::new(Vec::new(), 0.5, 0.5);
    let mut grn = grn_model::GrnModel::new_default();
    grn.load_from_file(path);

    grn.reset();
    grn.evolve(25);
    println!("{:?}", grn);

    // test model 
    // input concentrations 
    let mut y_res = Vec::<f64>::new();
    let mut x = Vec::<f64>::new();
    for step in 0..100
    {
        x.push((step as f64)/1.);
        grn.proteins[0].concentration = 0.125 * (step as f64)/100.;
        //grn.proteins[1].concentration =  0.125 * f64::sin(step as f64 + 1.0)/2.;
        grn.evolve(1);
        let y_hat = grn.proteins[1].concentration*2. -1.0 ;
        y_res.push(y_hat);
    }

    println!("res is {:?}", y_res);

    // do the plot things 

    let trace = Scatter::new(x, y_res).mode(Mode::Lines);

    let layout = Layout::new().x_axis(Axis::new().title(Title::from("X Axis")))
        .y_axis(Axis::new().title(Title::from("Y Axis")))
        .title(Title::from("My Plot"));

    let mut plot = Plot::new();
    plot.use_local_plotly();
    plot.add_trace(trace);
    plot.set_layout(layout);
    plot.show();
    //plot.write_html("out.html");

}

#[warn(dead_code)]
fn test_gene_id()
{
    // let rng : &dyn Rng;
    let first_gene = GrnGene::generate_random_gene(0, 0);
    let second_gene = GrnGene::generate_random_gene(0, 0);
    let third_gene: GrnGene = GrnGene::generate_random_gene(0, 0);
    let fourth_gene = GrnGene::generate_random_gene(0, 0);

    println!("fist gene id : {:?}", first_gene.get_id());
    println!("second gene id : {:?}", second_gene.get_id());
    println!("third_gene gene id : {:?}", third_gene.get_id());
    println!("fourth_gene gene id : {:?}", fourth_gene.get_id());
    

}
fn main() {
    test_grn();
}
