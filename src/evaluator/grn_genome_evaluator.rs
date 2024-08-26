use crate::evolver::{grn_gene::GrnGene, grn_genome::GrnGenome};


pub trait GrnGenomeEvaluator {

    fn evaluate(&self, gene : &GrnGenome) -> f64;
    // fn clone_box(&self) -> Box<dyn GrnGenomeEvaluator>;
}

pub trait CloneBox {
    fn clone_box(&self) -> Box<dyn GrnGenomeEvaluator>;
}

// Implement CloneBox for any type that implements your trait and Clone
impl<T> CloneBox for T
where
    T: 'static + Clone + GrnGenomeEvaluator, // Replace `Trait` with your actual trait name
{
    fn clone_box(&self) -> Box<dyn GrnGenomeEvaluator> {
        Box::new(self.clone())
    }
}
