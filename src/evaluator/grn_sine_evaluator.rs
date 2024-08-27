use super::grn_genome_evaluator::GrnGenomeEvaluator;
use crate::evolver::grn_genome::GrnGenome;

#[derive(Debug)]
pub struct SineEvaluator {
    num_grn_inputs : i32,
    num_grn_outputs : i32,
    num_evaluation : i32,
    name : String
}

impl SineEvaluator {
    fn new(input : i32, output : i32) -> Self
    {
        Self { num_grn_inputs: input, num_grn_outputs: output, num_evaluation: 0, name: "sine experience".to_string() }
    }
}

impl GrnGenomeEvaluator for SineEvaluator
{
    fn evaluate(&self, gene : &GrnGenome) -> f64 {
        gene.has_been_evaluated();
        0.0
    }
    

}

impl Clone for SineEvaluator {
    fn clone(&self) -> Self {
        Self { num_grn_inputs: self.num_grn_inputs, num_grn_outputs: self.num_grn_outputs, num_evaluation: 0, name: "sine experience".to_string() }
    }
}