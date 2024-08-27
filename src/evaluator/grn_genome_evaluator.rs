use std::string;

use crate::evolver::{grn_gene::GrnGene, grn_genome::GrnGenome};
use dyn_clone::DynClone;

pub trait GrnGenomeEvaluator : DynClone  {

    fn evaluate(&self, gene : &GrnGenome) -> f64;
    // fn clone_box(&self) -> Box<dyn GrnGenomeEvaluator>;
}

#[derive(Clone)]
pub struct defaultEvaluator {
    name : String,
}

impl GrnGenomeEvaluator for defaultEvaluator {
    fn evaluate(&self, gene : &GrnGenome) -> f64 {
        0.
    }
}


dyn_clone::clone_trait_object!(GrnGenomeEvaluator);

impl defaultEvaluator {
    pub fn new() -> Self
    {
        Self {name:"default".to_string()}
    }
}
